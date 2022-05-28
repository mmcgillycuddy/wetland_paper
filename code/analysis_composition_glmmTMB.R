library(glmmTMB)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(ochRe)
source("code/functions.R")
set.seed(1001)
###-----
###-----
### Composition models
###-----
###-----
# read in data
compo <- read.csv("data/derived/composition_data_long_March22.csv") %>% 
  mutate(Time_pt = as.numeric(Time_pt)) %>% 
  mutate(fire_treatment_aloc = fire_treatment) %>% 
  mutate(fire_treatment = ifelse(fire_treatment == "b" &  as.numeric(Time_pt) <= 6, "ub",fire_treatment)) %>% 
  mutate(fire_treatment = factor(fire_treatment, levels = c("ub","b"))) 

### Analysis dataset
### Remove sods and species with few presences
PA_data <- compo %>% 
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
### Traits model
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

### 
em = emmeans(comp_water_fire_fit, ~ Time_pt + water_treatment, at = list(Time_pt = c(1,9)))

base <- as.data.frame(diag(2*3))
combm = as.data.frame(em)[,1:2]
colnames(base)= rownames(base)=paste0(combm$water_treatment,combm$Time_pt)
selected_contrasts <- contrast(em, method =  list("H v. L by time interaction" =  ( (base$L9 - base$L1) - (base$H9 - base$H1)),
                                                  "H v. M by time interaction" =  ( (base$M9 - base$M1) - (base$H9 - base$H1))), 
                               adjust = "mvt", type = "response")
CI = confint(selected_contrasts) %>% as.data.frame()
1 - CI[,c(2,5,6)]

em = emmeans(comp_water_fire_fit, ~ fire_treatment + water_treatment, at = list(Time_pt = c(1,9)))
base <- as.data.frame(diag(2*3))
combm = as.data.frame(em)[,1:2]
colnames(base)= rownames(base)=paste0(combm$fire_treatment,combm$water_treatment)
selected_contrasts <- contrast(em, method = list("H v. L by b vs ub" =  (  (base$ubL - base$ubH) - (base$bL - base$bH)),
                                      "H v. M by b vs ub" = ( (base$ubM - base$ubH) - (base$bM - base$bH))), 
                               adjust = "mvt", type = "response")

CI = confint(selected_contrasts) %>% as.data.frame()
CI
### --------------------------------------
### Plot random effects - diag
### --------------------------------------
### Get 
object <- comp_water_fire_fit
cnms <- object$modelInfo$reTrms[["cond"]]$cnms   ## list of (named) terms and X columns
reStruc <- object$modelInfo$reStruc[[paste0("cond", "ReStruc")]] ## random-effects structure
flist <- object$modelInfo$reTrms[["cond"]]$flist ## list of grouping variables
levs <- lapply(flist, levels)

pl <- object$obj$env$parList(object$fit$par, object$fit$parfull)
nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
nb <- vapply(reStruc, function(x) x$blockReps, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)
bc <- vapply(reStruc, function(x) x$blockCode, numeric(1)) ## code block
nbseq <- rep.int(seq_along(nb), nb * nc)       ## splitting vector

block.plot <- 1 # diag re is the first random effect in model
bs <- ranef(object)
b <- as.matrix(bs$cond[[block.plot]])
b.plot.long <- as.data.frame(t(b)) %>%
  tibble::rownames_to_column("Coefficient") %>%
  tidyr::pivot_longer(cols = -Coefficient, names_to = names(bs$cond)[block.plot], values_to = "b" )
s1 <- TMB::sdreport(object$obj,  getJointPrecision = TRUE)
# subset first
s1.joint <- s1$jointPrecision
s1.b <- s1.joint[rownames(s1.joint)=="b", colnames(s1.joint)=="b"]
b.cols <- split(1:ncol(s1.b), nbseq)
b.cols.plot <-  b.cols[[block.plot]]
s1.b.plot <- s1.b[b.cols.plot, b.cols.plot]
s1.b.inv <- solve(s1.b.plot)
sd.b <- as.data.frame(sqrt(diag(s1.b.inv)))
coef.nam <- rep(cnms[[block.plot]], nb[[block.plot]])
nam.grp <- rep(levs[[block.plot]], each =nc[[block.plot]])
sd.re <- cbind(coef.nam, nam.grp, sd.b)
names(sd.re) <- c("Coefficient", names(levs[block.plot]), "sd.b")
ranef.plot <- left_join(b.plot.long, sd.re, by = intersect(names(b.plot.long), names(sd.re)))

load("results/comp_fire_water_plot_dataMarch22.RData")
### --------------------------------------
### Data for plots
### --------------------------------------
### ONly interested in interaction terms
highlight_coef = c("Time_pt:water_treatmentL", "water_treatmentL:fire_treatmentb",  "water_treatmentM:fire_treatmentb")
# add fixed effects
fixed_coefs <- fixef(comp_water_fire_fit)$cond %>% 
  as.data.frame() %>%
  rownames_to_column( var = "Coefficient") %>% 
  rename(fixed_coef = 2)
### 
ranef_fixed_plot <- ranef.plot %>%
  left_join(fixed_coefs, by = "Coefficient") %>% 
  mutate(conf.low = b - 1.96*sd.b,
         conf.high =  b + 1.96*sd.b,
         not.zero = ifelse(conf.high < 0 | conf.low > 0, 1, 0),
         comp_high = ifelse(Coefficient %in% highlight_coef &  not.zero, 1, 0 ),
         btot = b + fixed_coef,
         conf.low_tot = conf.low + fixed_coef,
         conf.high_tot = conf.high + fixed_coef,
         not.zero_tot = ifelse(conf.high_tot < 0 | conf.low_tot > 0, 1, 0),
         abund_high = ifelse(Coefficient %in% highlight_coef &  not.zero_tot, 1, 0 ),
         effect = ifelse(abund_high==1 & comp_high == 1,"abundance and composition",
                         ifelse(abund_high==1 & comp_high == 0,"abundance only",
                                ifelse(comp_high == 1, "composition only","neither"))) )

### plot.spp are the species that have at least one effect
plot.spp <- ranef_fixed_plot$Species[ranef_fixed_plot$effect!="neither"] 
# New facet label names for supp variable
supp.labs <- c("(a) Low against high \n water availability \n over time",
               "(b) Burnt against Unburnt sods \n in Low against High \n water availability",
               "(c) Burnt against Unburnt sods \n in Medium against High \n water availability")
names(supp.labs) <- highlight_coef

### Get full names for species 
sp_names <- read.csv("data/full_species_names_comp_analysis.csv") %>% 
  rename( Species_full = Full.species.name..in.italics.)
### --------------------------------------
### Treatment effect Plot
### Where treatment is fixed + random effect, (ie. beta + b = btot)
### --------------------------------------
### keep species that have an abundance effect
Spp_abundance_data <- ranef_fixed_plot %>% 
  filter(Species %in% plot.spp,
         Coefficient %in% highlight_coef,
         effect=="abundance only" | effect =="abundance and composition") %>%
  mutate(above_below = case_when(conf.high_tot < 0 ~ "below",
                                 conf.low_tot > 0 ~ "above")) %>% 
  droplevels() %>% 
  left_join(sp_names)

spp_abund_plot <- Spp_abundance_data %>% 
  ggplot(aes(x = Species_full, y = btot, color = above_below)) +
  geom_point()+ 
  geom_linerange(aes(ymin = conf.low_tot, ymax = conf.high_tot), size = 1, alpha = 0.5 ) +
  geom_hline(yintercept = 0) +
  facet_grid(~ Coefficient, scales = "free", 
             labeller = labeller(Coefficient = supp.labs), switch = "y") +
  coord_flip() + 
  theme_classic()  +
  scale_x_discrete(limits = rev) +
  xlab("Species with evidence of treatment effect") + 
  ylab("log odds ratio (strength of association)") +
  theme( axis.text = element_text( size = 12 ),
         axis.text.y = element_text( size = 10, face = "italic"),
         axis.text.x = element_text( size = 8),
         axis.title = element_text( size = 12),
         strip.text = element_text(size = 8),
         legend.position="none",
         strip.text.x = element_text(size=8, angle=0)) +
  scale_colour_manual(values = clrs2)


ggsave(file = "plots/spp_abund.png", spp_abund_plot, device = "png", height = 8, width = 8)


### --------------------------------------
### Composition  Plot
### i.e. diag(x | spp) random effect
### --------------------------------------
Spp_comp_data <- ranef_fixed_plot %>% 
  filter(Species %in% plot.spp,
         Coefficient %in% highlight_coef,
         effect == "composition only" | effect == "abundance and composition") %>% 
  mutate(above_below = case_when(conf.high_tot < fixed_coef ~ "below",
                                 conf.low_tot > fixed_coef ~ "above")) %>% 
  droplevels() %>% 
  left_join(sp_names)

spp_comp_plot <- Spp_comp_data %>% 
  ggplot(aes(x = Species_full, y = btot, color = above_below)) +
  geom_point()+ 
  geom_linerange(aes(ymin = conf.low_tot, ymax = conf.high_tot), size = 1, alpha = 0.5 ) +
  geom_hline(data = fixed_coefs %>% filter(Coefficient %in% highlight_coef), 
             aes(yintercept = fixed_coef), linetype  = "dashed", color = "gray60")+ 
  facet_grid(~ Coefficient, scales = "free", labeller = labeller(Coefficient = supp.labs), switch = "y") +
  coord_flip() + 
  theme_classic() +
  scale_x_discrete( limits = rev )+
  xlab("Species contributing to composition effect") + 
  ylab("log odds ratio (strength of association)") +
  theme( axis.text = element_text( size = 12 ),
         axis.text.y = element_text( size = 10, face = "italic"),
         axis.text.x = element_text( size = 8),
         axis.title = element_text( size = 12),
         strip.text = element_text(size = 8),
         legend.position="none",
         strip.text.x = element_text(size=8, angle=0)) +
  scale_color_manual(values = clrs2)

ggsave(file = "plots/spp_comb.png", spp_comp_plot, device = "png", height = 5, width = 7)

### --------------------------------------
### Traits
### Water trait
### --------------------------------------
pos = position_dodge(width=1)

em = emmeans(comp_water_trait_fit, ~ Time_pt*water_treatment*Water_trait_derived, at = list(Time_pt = c(1,9)),combine= TRUE)
base <- as.data.frame(diag(2*3*3))
combm = as.data.frame(em)[,1:3]

colnames(base)= rownames(base)=paste0("T",combm$Time_pt,combm$water_treatment,combm$Water_trait_derived)
selected_contrasts <- contrast(em, method = list("LH High " = ( (base$T9LHigh - base$T1LHigh) - (base$T9HHigh - base$T1HHigh)),  
                                                 "LH Medium" = ( (base$T9LMedium - base$T1LMedium) - (base$T9HMedium - base$T1HMedium)),
                                                 "LH Low" = ( (base$T9LLow - base$T1LLow)) - (base$T9HLow - base$T1HLow),
                                                 "MH High " = ( (base$T9MHigh - base$T1MHigh) - (base$T9HHigh - base$T1HHigh)),  
                                                 "MH Medium" = ( (base$T9MMedium - base$T1MMedium) - (base$T9HMedium - base$T1HMedium)),
                                                 "MH Low" = ( (base$T9MLow - base$T1MLow)) - (base$T9HLow - base$T1HLow)), 
                               adjust = "mvt")


loc_labs  =c("(a) Low against high \n water availability \n over time", "(a) Medium against high \n water availability \n over time") 

water_trait_mm <- selected_contrasts %>% 
  as.data.frame() %>% 
  mutate(contr = substr(contrast,1,2),
         Trait = str_trim(substr(contrast,4,9))) %>% 
  mutate(water_level = factor(contr, levels = c( "LH", "MH"), labels = loc_labs),
         Trait = factor(Trait, levels = c("High","Medium","Low"), labels = c( "Strong hydrophile",  "Moderate hydrophile", "Non hydrophile" ))) %>% 
  ggplot( aes(x = estimate, y = Trait)) + 
  geom_point()+
  geom_segment(aes(x = estimate - 1.96*SE, xend = estimate + 1.96*SE, 
                   y = Trait, yend = Trait, colour = Trait, size = 1,
                   alpha = 0.7))+
  theme_classic()+
  geom_vline(xintercept = 0)+
  facet_grid(~ water_level)+
  theme(legend.position = "none")+
  xlab("log odds ratio (strength of association)") +
  ylab("Water requirement") +
  scale_color_manual(values = clrs3) 

ggsave(file = "plots/water_trait.png", water_trait_mm, device = "png", height = 4, width = 8)
### --------------------------------------
### Fire trait
### --------------------------------------
pos = position_dodge(width=1)
em_fire = emmeans(comp_fire_trait_fit, ~ fire_treatment*water_treatment*Fire.response, data = PA_data, combine= TRUE)
base <- as.data.frame(diag(2*3*3))
combm = as.data.frame(em_fire)[,1:3]

colnames(base)= rownames(base) = paste0(combm$fire_treatment,combm$water_treatment,combm$Fire.response)
select_cont_fire <- contrast(em_fire, method = list("LH Killed" =  ( (base$bLKilled - base$ubLKilled) - (base$bHKilled - base$ubHKilled)),  
                                                 "LH Resprouter" =  (  (base$bLResprouter - base$ubLResprouter) - (base$bHResprouter - base$ubHResprouter)),
                                                 "LH Other" =  ( (base$bLOther  - base$ubLOther ) - (base$bHOther  - base$ubHOther )),
                                                 "MH Killed" = ( (base$bMKilled - base$ubMKilled) - (base$bHKilled - base$ubHKilled)),  
                                                 "MH Resprouter" =  ( (base$bMResprouter - base$ubMResprouter) - (base$bHResprouter - base$ubHResprouter)),
                                                 "MH Other" = ( (base$bMOther  - base$ubMOther ) - (base$bHOther  - base$ubHOther ))),
                               adjust = "mvt")
fire_trait_mm <- select_cont_fire %>% 
  as.data.frame() %>% 
  mutate(contr = substr(contrast,1,2),
         Trait = str_trim(substr(contrast,4,13))) %>% 
  mutate(water_level = factor(contr, levels = c( "LH", "MH"), labels = supp.labs[-1]),
         Trait = factor(Trait, levels = c("Killed","Resprouter","Other"))) %>% 
  ggplot( aes(x = estimate, y = Trait)) + 
  geom_point() +
  geom_segment(aes(x = estimate - 1.96*SE, xend = estimate + 1.96*SE, 
                   y = Trait, yend = Trait, colour = Trait, size = 1, alpha = 0.7)) +
  theme_classic() +
  geom_vline(xintercept = 0) +
  facet_grid( ~ water_level) +
  theme(legend.position = "none") +
  xlab("log odds ratio (strength of association)") +
  ylab("Fire response") +
  scale_color_manual(values = clrs3) 

ggsave(file = "plots/fire_trait.png", fire_trait_mm, device = "png", height = 4, width = 8)

save.image("results/comp_fire_water_plot_dataMarch22.RData")
