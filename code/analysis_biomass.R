library(dplyr)
library(ggplot2)
library(tidyr)
library(lmerTest)
library(lme4)
library(emmeans)
library(kableExtra)
source("code/functions.R")
### --------------------------------------
### Data and plotting of biomass
### --------------------------------------
times <- read.csv("data/time.csv") %>% 
  select(-Time_pt) %>% 
  mutate(Median_time_factor = factor(Median_time)) 

baci_bio <- read.csv("data/baci_bio_df3_sod_live_altered.csv") %>% 
  mutate(Sod_biomass=as.numeric(Sod_biomass)) %>% 
  mutate(water_treatment = factor(water_treatment, levels = c("H", "M", "L"))) %>% 
  mutate(fire_treatment_aloc = fire_treatment) %>% 
  mutate(fire_treatment = ifelse(fire_treatment == "b" &  Median_time_pt == 592, "ub",fire_treatment)) %>% 
  mutate(fire_treatment = factor(fire_treatment, levels = c("ub","b"))) %>% 
  left_join(times, by = "Median_time_pt") %>%  
  mutate(Median_time_factor = factor(Median_time)) 
## need to add other time points to the data
blanks <- expand.grid(water_treatment = levels(baci_bio$water_treatment),
                      fire_treatment = levels(baci_bio$fire_treatment),
                      Median_time = unique(times$Median_time),
                      veg = unique(baci_bio$veg)) 
df_to_plot <- bind_rows(blanks,baci_bio) 

## Box plot of raw data
p <- ggplot(df_to_plot, aes(x = Median_time, y = (Sod_biomass + 1), fill = water_treatment, color = water_treatment, group = interaction(water_treatment, Median_time))) + 
  geom_boxplot() +
  theme_classic()+
  facet_grid( veg ~ fire_treatment)+
  theme(legend.position = "bottom")+
  xlim(0, 1500)+
  scale_y_continuous(breaks = c(2,11,101), labels = c(1,10,100), trans = "log10")

### --------------------------------------
### Model: Linear mixed model
### --------------------------------------
bio_live <- lmer(log(Sod_biomass+1) ~ Swamp + veg + Median_time_factor*water_treatment +  
                   water_treatment*fire_treatment + Median_time_factor*fire_treatment + (1 | Sod), data = baci_bio)
### residuals look good
plot(residuals(bio_live) ~ fitted(bio_live, re.form = NULL))
### Type III anova test 
anova(bio_live, type = "III")

###  Plot of water_treatment over time effect in burnt and unburnt sods
emmeans_plot = emmip(bio_live,  ~ water_treatment ~ Median_time_factor| fire_treatment , CIs = TRUE, type = "response") 
joined <- emmeans_plot$data %>% 
  left_join(times, by = "Median_time_factor")
## need to add other time points to the data
blanks <- expand.grid(water_treatment = levels(joined$water_treatment),
                      fire_treatment = levels(joined$fire_treatment),
                      Median_time = unique(times$Median_time))
plot_df <- joined %>% 
  bind_rows(blanks)
pos = position_dodge(width = 80)

bio_plot <- plot_df %>%
  ggplot(aes(Median_time, yvar, color = water_treatment, shape = fire_treatment, linetype = fire_treatment), position = position_dodge2(width = 0.1)) +
  geom_point( position = pos, size = 2) +
  geom_path(position = pos, size = 1) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), position = pos,  alpha = 0.5, width = 0, size = 2) +
  theme_classic() +
  xlim(0,1500) +
  xlab("Time since experiment commenced (days)") + 
  ylab("Mean biomass (+/- 95% CI) per mesocosm (g)") +
  theme(legend.position = "none",
        axis.text = element_text( size = 12), 
        axis.title = element_text( size = 14)) +
  scale_colour_manual(values = clrs3)
bio_plot
ggsave(plot = bio_plot, file = "plots/biomass_plot.tiff" , width = 150, height = 150, units = "mm", device = "tiff")

### --------------------------------------
# Post hoc analysis of time * water interaction unburnt sods
# Look at pairwise comparisons of the water treatment differences change over time
### --------------------------------------
em_water <- emmeans(bio_live,  ~ water_treatment + Median_time_factor + fire_treatment, combine= TRUE )
base <- as.data.frame(diag(3*2*2))
combm = as.data.frame(em_water)[,1:3]

colnames(base)= rownames(base)=paste0(combm$water_treatment,combm$Median_time_factor,combm$fire_treatment)
selected_contrasts <- contrast(em_water, method =  list("H v. M by time interaction" = (base$H1261ub - base$M1261ub) - (base$H587ub - base$M587ub),
                                      "H v. L by time interaction" =  (base$H1261ub - base$L1261ub) - (base$H587ub - base$L587ub),
                                      "M v. L by time interaction" =  (base$M1261ub - base$L1261ub) - (base$M587ub - base$L587ub)), 
                               adjust = "mvt", type = "response")

## Condifidence intervals of contrasts
cbind(summary(selected_contrasts), confint(selected_contrasts)[, 5:6]) %>% 
  as.data.frame() %>% 
  mutate_at(2:8, round, 3)
