library(tidyr)
library(dplyr)

treatment <- read.csv("data/water_fire_treatment.csv") %>% 
  select(Sod,water_treatment, fire_treatment, veg)

species_data <- read.csv("data/long_comp_2017_2020_Mar_2022_update.csv", strip.white = TRUE) %>% 
  select(Sod, Species, Time_pt, HvMax, live_dead) %>% 
  mutate(HvMax = as.numeric(HvMax))
# filter out the dead observations ("d"), "" and "dying" is included
species_data <- species_data %>% 
  filter(live_dead != "d") %>% 
  select(-live_dead)

species_names <- read.csv("data/species_list_glasshouse_updated_March_22.csv", header = FALSE)
# start with a line for each Species by time by Sod combination and all the species data
expanded_data<- expand.grid(Species = species_names$V1,  Time_pt = 1:9, Sod = sort(unique(species_data$Sod))) %>% 
  full_join(species_data) %>% 
  filter(Species %in% species_names$V1) 

compo = expanded_data %>% 
  mutate(PA_each = ifelse(!is.na(HvMax),1,0)) %>% 
  group_by(Sod, Time_pt, Species) %>% 
  mutate(PA = ifelse(any(PA_each == 1), 1, 0)) %>% 
  distinct(Sod, Time_pt, Species, .keep_all = TRUE) %>% 
  select(-HvMax, - PA_each) %>% 
  mutate(Swamp = substr(Sod, 1,2)) %>% 
  left_join(treatment) 

write.csv(compo,"data/derived/composition_data_long_March22.csv")

### Make data wide
# a. Change burnt to unburnt if it is before time point 6
compo_tmp <- compo %>% 
  mutate(Time_pt = as.factor(Time_pt)) %>% 
  mutate(fire_treatment_aloc = fire_treatment) %>% 
  mutate(fire_treatment = ifelse(fire_treatment == "b" &  as.numeric(Time_pt) <= 6, "ub",fire_treatment)) %>% 
  mutate(fire_treatment = factor(fire_treatment, levels = c("ub","b"))) 

compo_wide <- compo_tmp %>% 
    pivot_wider(names_from = Species, values_from = PA, id_cols = c( Time_pt, Sod, Swamp:fire_treatment_aloc) )

write.csv(compo_wide,"data/derived/composition_data_wide_March22.csv")


### Richness dataset
times <- read.csv("data/time.csv") %>% 
  select(-Median_time_pt) %>% 
  mutate(Time_pt_factor = factor(Time_pt))

treatment <- read.csv("data/water_fire_treatment.csv") %>% 
  select(Sod,water_treatment, fire_treatment, veg)

rich_data <- read.csv("data/derived/composition_data_long_March22.csv") %>% 
  group_by(Sod, Time_pt) %>% 
  summarise(richness = sum(PA), Swamp = first(Swamp)) %>% 
  left_join(treatment) %>% 
  mutate(fire_treatment_aloc = fire_treatment) %>% 
  mutate(fire_treatment = ifelse(fire_treatment == "b" &  Time_pt  <= 6, "ub",fire_treatment)) %>% 
  mutate(water_treatment = factor(water_treatment, levels = c("H", "M", "L"))) %>% 
  mutate(fire_treatment = factor(fire_treatment, levels = c("ub","b"))) %>% 
  mutate(Time_num = Time_pt) %>% 
  mutate(Time_pt_factor = factor(Time_pt)) %>% 
  left_join(times)

write.csv(rich_data,"data/derived/richness_long.csv")

