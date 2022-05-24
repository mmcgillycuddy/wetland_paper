library(glmmTMB)
library(tidyverse)
library(ggplot2)
source("functions.R")
set.seed(1001)
###-----
###-----
### Composition models
###-----
###-----
# read in data
compo <- read.csv("data/derived/composition_data_long_March22.csv")

# Change burnt to unburnt if it is before time point 6
PA_long <- compo %>% 
  mutate(Time_pt = as.numeric(Time_pt)) %>% 
  mutate(fire_treatment_aloc = fire_treatment) %>% 
  mutate(fire_treatment = ifelse(fire_treatment == "b" &  as.numeric(Time_pt) <= 6, "ub",fire_treatment)) %>% 
  mutate(fire_treatment = factor(fire_treatment, levels = c("ub","b"))) 

### Analysis dataset
### Remove sods and species with few presences
PA_data <- PA_long %>% 
  group_by(Species) %>% 
  mutate(n_pres_spp = sum(PA)) %>% 
  ungroup(Species) %>% 
  group_by(Sod) %>% 
  mutate(n_pres_sod = sum(PA)) %>% 
  ungroup(Sod) %>% 
  filter(n_pres_sod > 9 & n_pres_spp > 9) %>% 
  droplevels() 

## Traits data
trait_data <- read.csv("data/water_fire_response_species.csv") %>% 
  select(Species.code, Water_trait_derived, Fire.response  ) %>% 
  rename(Species = Species.code)

PA_trait_data <- left_join(PA_data, trait_data, by = "Species")
PA_data <- PA_trait_data
table(PA_data$Species, PA_data$PA)
table(PA_data$Sod, PA_data$PA)
nspp <- nlevels(as.factor(PA_data$Species))
nsods <- nlevels(as.factor(PA_data$Sod))
###-----
### Composition testing, i.e., testing random effect
### Variables of interest: 
### - a) water_treatment*fire_treatment
### - b) Time_pt*water_treatment
###-----
###-----
### Null Model for 1a: comp_water.fire_fixed_fit
### water_treatment*fire_treatment is fixed, i.e., not in random effect
###-----
fmla_water.fire_fixed = as.formula(PA ~ Swamp + veg +  Time_pt*water_treatment +  Time_pt*fire_treatment + water_treatment*fire_treatment + 
                                     diag(Swamp + veg +  Time_pt*water_treatment + Time_pt*fire_treatment | Species) + 
                                     rr(0 + Species | Sod, 2))

time.start <- Sys.time()
comp_water.fire_fixed_fit <- glmmTMB(fmla_water.fire_fixed,
                                     family = binomial(), 
                                     control = glmmTMBControl(start_method = list(method = "res")),
                                     data = PA_data)
time.end <- Sys.time()
time_water.fire_fixed <- time.end - time.start 
###-----
### Null Model for 1b: comp_time.water_fixed_fit
###-----
fmla_time.water_fixed = as.formula(PA ~ Swamp + veg +  Time_pt*water_treatment + Time_pt*fire_treatment + water_treatment*fire_treatment + 
                                     diag(Swamp + veg +  Time_pt*fire_treatment + water_treatment*fire_treatment | Species) + 
                                     rr(0 + Species | Sod, 2))
### Get the structure of the model
comp_time.water_fixed_nofit <- glmmTMB(fmla_time.water_fixed, family = binomial(), doFit = F, data = PA_data)
# Get initial parameter estimates from the fitted model comp_water.fire_fixed_fit
start.par.m1b <- startParFitted(comp_water.fire_fixed_fit, comp_time.water_fixed_nofit)

time.start <- Sys.time()
comp_time.water_fixed_fit <- glmmTMB(fmla_time.water_fixed, family = binomial(), data = PA_data,
                               start = list(beta = start.par.m1b$start.beta, theta = start.par.m1b$start.theta, b = start.par1b$start.b))
time.end <- Sys.time()
time_time.water_fixed <- time.end - time.start 

###-----
### Model 2: Full model with all effects (i.e., variables of interest in both fixed and random effects)
###-----
fmla_water_fire = as.formula(PA ~ Swamp + veg +  Time_pt*water_treatment + Time_pt*fire_treatment + water_treatment*fire_treatment + 
                               diag(Swamp + veg +  Time_pt*water_treatment + Time_pt*fire_treatment + water_treatment*fire_treatment | Species) + 
                               rr(0 + Species | Sod, 2))
### Initialise starting parameters
### Getting the structure of the parameters from the model without fitting it
comp_water_fire_nofit <- glmmTMB(fmla_water_fire, family = binomial(), data = PA_data, doFit = F)
start.par.m2 <- startParFitted(comp_water.fire_fixed_fit , comp_water_fire_nofit)

### Full model fit: comp_water_fire_fit 
time.start.water_fire <- Sys.time()
comp_water_fire_fit <- glmmTMB(fmla_water_fire, family = binomial(), data = PA_data, 
                               start = list(beta = start.par.m2$start.beta, theta = start.par.m2$start.theta, b = start.par.m2$start.b))
time.end.water_fire <- Sys.time()
time.water_fire <- time.end.water_fire - time.start.water_fire 
###-----
### test random effects
###-----
anova(comp_water_fire_fit, comp_water.fire_fixed_fit)
anova(comp_water_fire_fit, comp_time.water_fixed_fit)
###-----
###-----
### Traits model
###-----
##-----
## Water trait
##-----
fmla_water_trait = as.formula(PA ~ Swamp + veg +  Time_pt*water_treatment + Time_pt*fire_treatment + water_treatment*fire_treatment + 
                                Time_pt*water_treatment*Water_trait_derived +
                                diag(Swamp + veg +  Time_pt*water_treatment + Time_pt*fire_treatment + water_treatment*fire_treatment | Species) + 
                                rr(0 + Species | Sod, 2))
### Initialise starting parameters
###-----
### Getting the structure of the parameters from the model without fitting it
comp_trait_nofit <- glmmTMB(fmla_water_trait, family = binomial(), data = PA_data, doFit = F)
# Get initial parameter estimates from the fitted model comp_water.fire_fixed_fit
start.par.wt <- startParFitted(comp_water.fire_fixed_fit, comp_trait_nofit)

### Trait water model fit
comp_water_trait_fit <- glmmTMB(fmla_water_trait, family = binomial(), data = PA_data, 
                                start = list(beta = start.par.wt$start.beta, theta = start.par.wt$start.theta, b = start.par.wt$start.b))


##-----
## Fire trait
##-----
fmla_fire_trait = as.formula(PA ~ Swamp + veg +  Time_pt*water_treatment + Time_pt*fire_treatment + water_treatment*fire_treatment + 
                                water_treatment*fire_treatment*Fire.response + 
                                diag(Swamp + veg +  Time_pt*water_treatment + Time_pt*fire_treatment + water_treatment*fire_treatment | Species) + 
                                rr(0 + Species | Sod, 2))
### Initialise starting parameters
comp_fire_trait_nofit <- glmmTMB(fmla_fire_trait, family = binomial(), data = PA_data, doFit = F)
start.par.ft <- startParFitted(comp_water.fire_fixed_fit, comp_fire_trait_nofit)

### Fire trait model fit
comp_fire_trait_fit <- glmmTMB(fmla_fire_trait,
                                family = binomial(), 
                                start = list(beta = start.beta, theta = start.theta, b = start.b),
                                data = PA_data)

## Test traits
anova(comp_water_fire_fit, comp_water_trait_fit)
anova(comp_water_fire_fit, comp_fire_trait_fit)

save(comp_fire_trait_fit, file = "results/comp_fire_trait_model.RData")
load( "results/comp_fire_water_plot_dataMarch22.RData")
