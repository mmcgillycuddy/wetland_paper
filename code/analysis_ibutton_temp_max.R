library(dplyr)
library(ggplot2)
library(tidyverse)
library(lmerTest)
library(lme4)
library(emmeans)
library(glmmTMB)
source("code/functions.R")
ibutton_temps <- read.csv("data/ibutton_max_temps_1.csv") %>% 
  mutate(Soil_depth = factor(Soil_depth))

## Model and residual checks
ibutton_max.lm <- glmmTMB(log(Max_temp) ~ Water_level*Soil_depth +(1|Sod), data = ibutton_temps)
scatter.smooth(residuals(ibutton_max.lm)~fitted(ibutton_max.lm, re.form = NULL))
scatter.smooth(sqrt(abs(residuals(ibutton_max.lm)))~fitted(ibutton_max.lm, re.form = NULL), col = "red")

## Hypothesis test
car::Anova(ibutton_max.lm, type = "II", test.statistic = "Chisq")
## Plot of estimated model

temp_plot <- emmip(ibutton_max.lm,  ~ Water_level ~ Soil_depth , CIs = TRUE, type = "response")+
  xlab("Soil depth (cm)")+
  ylab("Mean maximum temperature (°C) \n (+/- 95% CI)")+
  ylim(0,100)+
  theme(axis.text = element_text(size = 12))+
  theme(text=element_text(size=14))+
  scale_color_manual(values = c(clrs3[1], clrs3[3], clrs3[2]))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none") 
temp_plot

ggsave(plot = temp_plot, file = "plots/temp_plot.tiff" , width = 110, height = 100, units = "mm", device = "tiff")

depth_post <- emmeans(ibutton_max.lm,  ~  Water_level +  Soil_depth, combine= TRUE )
base <- as.data.frame(diag(6))
combm = as.data.frame(depth_post)[,1:4]

colnames(base)= rownames(base)=paste0(combm$Water_level,combm$Soil_depth)

selected_contrasts <- contrast(depth_post, method = 
                                 list("L1/M1" = (base$L1 - base$M1),
                                      "L1/H1" = (base$L1 - base$H1),
                                      "M1/H1" = (base$M1 - base$H1),
                                      "L3/M3" = (base$L3 - base$M3),
                                      "L3/H3" = (base$L3 - base$H3),
                                      "M3/H3" = (base$M3 - base$H3)),
                               adjust = "mvt", type = "response")

cbind(summary(selected_contrasts),confint(selected_contrasts)[,5:6]) %>% 
  as.data.frame() %>% 
  mutate_at(2:8, round, 3)

plot(selected_contrasts, xlab = "max temp diff")+
  geom_vline(xintercept = 1)+
  theme_classic()


