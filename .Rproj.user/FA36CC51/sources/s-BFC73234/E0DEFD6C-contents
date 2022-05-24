library(glmmTMB)
library(tidyr)
library(dplyr)
library(ggplot2)
library(statmod)
library(ggpubr)

# read in data

treatment <- read.csv("data/water_fire_treatment.csv") %>% 
  select(Sod, water_treatment, fire_treatment, veg) %>% 
  mutate(Swamp = substr(Sod,1,2))

compo <- read.csv("data/biomass_by_species_wide_2_updated.csv") %>% 
  select(-X) %>% 
  left_join(treatment) %>% 
  filter(!(fire_treatment == "ub" & Time_pt == 6)) %>% 
  mutate(fire_treatment = ifelse(fire_treatment == "b" &  Time_pt == 6, "ub",fire_treatment)) 

### Fix species columns that were recorded twice with different spellings
### remove one column full of zeros
all(compo$Xyris.oper==0)
compo <- select(compo, -Xyris.oper)
### combine Lepy.scar and lepi.scar
compo <- compo %>% 
  mutate(Lepi.scar = Lepi.scar, 
         Lepi.scar = ifelse(Lepy.scar!=0 & Lepi.scar==0, Lepy.scar, Lepi.scar) ) %>% 
  select(-Lepy.scar)

compo_long<- compo %>% 
  pivot_longer(Bank.eric:Acac.ptyc, names_to = "Species", values_to = "biomass") %>% 
  mutate(PA = biomass>0)

## Water only, time point 6 only
## is there a difference between water treatments at time 6, before the fire

sods_low <- compo_long %>% 
  filter(Time_pt == 6) %>% 
  group_by(Sod) %>% 
  summarise(rich = sum(PA)) %>% 
  filter(rich <= 1)

species_low <- compo_long %>% 
  filter(Time_pt == 6) %>% 
  group_by(Species) %>% 
  summarise(rich = sum(PA)) %>% 
  filter(rich < 10)

compo_mod <- compo_long %>% 
  filter(Time_pt == 6) %>% 
  group_by(Species) %>% 
  mutate(n_pres_spp = sum(PA)) %>% 
  ungroup(Species) %>% 
  group_by(Sod) %>% 
  mutate(n_pres_sod = sum(PA)) %>% 
  ungroup(Sod) %>% 
  #remove sods and species with few presences
  filter(n_pres_sod > 9 & n_pres_spp > 9) %>% 
  droplevels() 

#for time * water
formula_base = as.formula(biomass ~ Swamp + water_treatment + diag(Swamp|Species) + rr(Species + 0|Sod,2))
formula_time_water = as.formula( biomass ~ Swamp + water_treatment + diag(Swamp  + water_treatment|Species) + rr(Species + 0|Sod, 2))

st <- Sys.time()
mod_base <- glmmTMB(formula_base, data = compo_mod,
                    control = glmmTMBControl(start_method = list(method = "res")),
                       family=tweedie())

mod_time_water <-glmmTMB(formula_time_water, data = compo_mod,
                         control = glmmTMBControl(start_method = list(method = "res")),
                         family=tweedie())
Sys.time() - st

time_water_aov <- anova(mod_base, mod_time_water)#
time_water_aov

library(DHARMa)
mod_base_res <- simulateResiduals(mod_base)
mod_water_res <- simulateResiduals(mod_time_water)


# save(time_water_aov,file = "results/water_6")
load(file = "results/water_6")

## Water fire interaction, time point 9 only
sods_low <- compo_long %>% 
  filter(Time_pt == 9) %>% 
  group_by(Sod) %>% 
  summarise(rich = sum(PA)) %>% 
  filter(rich <= 1)

species_low <- compo_long %>% 
  filter(Time_pt == 9) %>% 
  group_by(Species) %>% 
  summarise(rich = sum(PA)) %>% 
  filter(rich < 7)
compo_mod <- compo %>% 
  select(-all_of(species_low$Species)) %>% 
  filter(!Sod %in% sods_low$Sod)

#for time * water
mv_formula_noint = Y   ~  Swamp + veg + water_treatment + fire_treatment

mv_formula_full = Y   ~  Swamp + veg + water_treatment*fire_treatment

st <- Sys.time() 

mv_mod_noint <- manyany("glm", Y, mv_formula_noint, data=X,
                        family=statmod::tweedie(var.power=1.1, link.power=0), var.power=1.1)

mv_mod_full <- manyany("glm", Y, mv_formula_full, data=X,
                       family=statmod::tweedie(var.power=1.1, link.power=0), var.power=1.1)



Sys.time() - st



water_fire_aov <- anova(mv_mod_noint, mv_mod_full, nBoot = 999, p.uni	= "unadjusted")#


save(water_fire_aov,file = "results/water_fire")

load(file = "results/water_fire")
water_fire_aov




