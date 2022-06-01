rm(list=ls())
library(dplyr)
library(plyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library (plotly)
library(lme4)
library(emmeans)
library(lmerTest)
source("functions.R")
treatment <- read.csv("data/water_fire_treatment_1.csv") %>% 
  select(Sod, Swamp, water_treatment, fire_treatment, veg) 

times <- read.csv("data/time_moisture.csv") %>% 
  select(-Median_time_pt) %>% 
  mutate(Time_pt_factor = factor(Time_pt))
sod_moisture<-read.csv("data/sod_moisture_all_dates.csv", header=TRUE) %>% 
  mutate(moist_1=as.numeric(moist_1)) %>% 
  mutate(moist_2=as.numeric(moist_2)) %>% 
  mutate(moist_3=as.numeric(moist_3)) %>% 
  mutate(Date = dmy(Date))

dates = sod_moisture %>% 
  group_by(Time_pt) %>% 
  summarise(midtime = median(Date,na.rm = TRUE))

sod_moisture = left_join(sod_moisture, dates, by = "Time_pt")%>% 
  rowwise() %>% 
  mutate(moist_av = mean(c(moist_1, moist_2, moist_3), na.rm = T))

moist <- sod_moisture %>% 
  group_by(Sod,Time_pt) %>% 
  left_join(treatment) %>% 
  mutate(fire_treatment_aloc = fire_treatment) %>% 
  mutate(fire_treatment = ifelse(fire_treatment == "b" &  Time_pt  <= 11, "ub",fire_treatment)) %>% ### check - time 1-11 are unburnt and 12 onwards are burnt
  mutate(water_treatment = factor(water_treatment, levels = c("H", "M", "L"))) %>% 
  mutate(fire_treatment = factor(fire_treatment, levels = c("ub","b"))) %>% 
  mutate(Time_num = Time_pt) %>% 
  mutate(Time_pt_factor = factor(Time_pt)) %>% 
  left_join(times)

pos = position_dodge(width=60)

p <- ggplot(moist, aes(x = Median_time, y=moist_av, fill = water_treatment, color = water_treatment)) + 
  # geom_boxplot(aes(group = interaction(water_treatment, Median_time)), size = 0.8, alpha = 0.4) +
  theme_classic()+
  geom_smooth(se= FALSE)+
  geom_jitter(position=pos, size = 0.8, alpha = 0.4)+
  ylab("soil moisture")+
  facet_grid(veg~fire_treatment)+
  theme(legend.position = "bottom")
p

### Model
#### Fit model and check assumptions
moist_time <- lmer(moist_av ~ Swamp + veg +  Time_pt_factor*water_treatment +  water_treatment*fire_treatment + Time_pt_factor*fire_treatment + (1 | Sod), data = moist)

scatter.smooth(residuals(moist_time)~fitted(moist_time, re.form = NULL))
scatter.smooth(sqrt(abs(residuals(moist_time)))~fitted(moist_time, re.form = NULL), col = "red")

### Inference

anova(moist_time, type = "III")

### Post hoc
#### Post hoc plot. 
emmeans_plot = emmip(moist_time,  ~ water_treatment ~ Time_pt_factor| fire_treatment , CIs = TRUE)+
  theme_classic()

## plot
emmeans_plot
joined <- emmeans_plot$data %>%
  left_join(times, by = "Time_pt_factor")
blanks <- expand.grid(water_treatment = levels(joined$water_treatment),
                      fire_treatment = levels(joined$fire_treatment),
                      Median_time = unique(times$Median_time))

plot_df <- joined %>% 
  bind_rows(blanks)

pos = position_dodge(width=70)
plot_moist <- plot_df %>% 
  ggplot(aes(Median_time,yvar, color = water_treatment, shape = fire_treatment, linetype = fire_treatment))+
  geom_point( position=pos)+
  geom_path(position=pos)+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), position=pos, alpha = 0.5, width = 0, lwd = 1.6)+
  xlab("Time since experiment commenced (days)")+
  ylab("Mean % volume soil moisture (+/- 95% CI)")+
  guides(fill=guide_legend(title=NULL))+
  xlim(0,1500)+
  ylim(-5,80)+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))+
  theme(plot.title = element_text(size = 14))+
  theme(strip.text = element_text(size = 14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = clrs3)

ggsave(plot = plot_moist, file = "plots/moisture_plot.tiff" , width = 180, height = 150, units = "mm", device = "tiff")

### burning increases soil moisture
em_fire_effect <- emmeans(moist_time,  ~ fire_treatment, 
                          at =  list(Time_pt_factor = c("12","13","14","15","16")), combine= TRUE)
p = pairs(em_fire_effect)
confint(p)


