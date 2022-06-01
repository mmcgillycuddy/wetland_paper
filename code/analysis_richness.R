library(dplyr)
library(ggplot2)
library(tidyr)
library(lmerTest)
library(lme4)
library(emmeans)
library(Cairo)
library(kableExtra)
library(ochRe)

treatment <- read.csv("data/water_fire_treatment.csv") %>% 
  select(Sod,water_treatment, fire_treatment, veg) 

times <- read.csv("data/time.csv") %>% 
  select(-Median_time_pt) %>% 
  mutate(Time_pt_factor = factor(Time_pt))

rich_data <- read.csv("data/derived/composition_data_long_March22.csv") %>% 
  group_by(Sod,Time_pt) %>% 
  summarise(richness = sum(PA), Swamp = first(Swamp)) %>% 
  left_join(treatment) %>% 
  mutate(fire_treatment_aloc = fire_treatment) %>% 
  mutate(fire_treatment = ifelse(fire_treatment == "b" &  Time_pt  <= 6, "ub",fire_treatment)) %>% 
  mutate(water_treatment = factor(water_treatment, levels = c("H", "M", "L"))) %>% 
  mutate(fire_treatment = factor(fire_treatment, levels = c("ub","b"))) %>% 
  mutate(Time_num = Time_pt) %>% 
  mutate(Time_pt_factor = factor(Time_pt)) %>% 
  left_join(times)

pos = position_dodge(width=60)

p <- ggplot(rich_data, aes(x = Median_time, y = richness, fill = water_treatment, color = water_treatment)) + 
  # geom_boxplot(aes(group = interaction(water_treatment, Median_time)), size = 0.8, alpha = 0.4) +
  theme_classic()+
  geom_smooth(se= FALSE)+
  geom_jitter(position=pos, size = 0.8, alpha = 0.4)+
  ylab("Richness")+
  facet_grid(veg~fire_treatment)+
  theme(legend.position = "bottom")
### --------------------------------------
### Model: Linear mixed model
### --------------------------------------
rich_live <- lmer(richness ~ Swamp + veg +  Time_pt_factor*water_treatment + Time_pt_factor*fire_treatment + water_treatment*fire_treatment  + (1 | Sod),
                  data = rich_data)
plot(residuals(rich_live)~fitted(rich_live, re.form = NULL))
anova(rich_live, type = "III")

## Post hoc plots
emmeans_plot = emmip(rich_live,  ~ water_treatment ~ Time_pt_factor| fire_treatment ,
                     at = list(Time_pt_factor = unique(rich_data$Time_pt_factor)),
                     CIs = TRUE) 

joined <- emmeans_plot$data %>%
  mutate(Time_pt_factor = as.factor(Time_pt_factor)) %>% 
  left_join(times, by = "Time_pt_factor")
blanks <- expand.grid(water_treatment = levels(joined$water_treatment),
                      fire_treatment = levels(joined$fire_treatment),
                      Median_time = unique(times$Median_time))
plot_df <- joined %>% 
  bind_rows(blanks)
pos = position_dodge(width=70)

rich_plot <- plot_df %>% 
  ggplot(aes(Median_time, yvar, color = water_treatment, shape = fire_treatment, linetype = fire_treatment))+
  geom_point( position = pos, size = 2 ) +
  geom_path( position = pos , size = 1) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), position = pos, alpha = 0.5, width = 0, size = 2) +
  theme_classic() +
  xlab("Time since experiment commenced (days)") + 
  ylab("Mean species richness (+/- 95% CI) per mesocosms (g)") +
  theme(legend.position = "none",
        axis.text = element_text( size = 12), 
        axis.title = element_text( size = 14)) +
  scale_colour_manual(values = clrs3)

ggsave(plot = rich_plot, file = "plots/richness_plot.tiff" , width = 180, height = 150, units = "mm", device = "tiff")

### --------------------------------------
### Post hoc analysis of interactions
### --------------------------------------
em_water <- emmeans(rich_live,  ~ water_treatment + Time_pt_factor + fire_treatment, combine= TRUE )
base <- as.data.frame(diag(3*9*2))
combm = as.data.frame(em_water)[,1:3]

colnames(base) = rownames(base) = paste0(combm$water_treatment,combm$Time_pt_factor,combm$fire_treatment)
selected_contrasts <- contrast(em_water, method = 
                                 list("H v. M by time interaction" = (base$H9ub - base$M9ub) - (base$H6ub - base$M6ub),
                                      "H v. L by time interaction" = (base$H9ub - base$L9ub) - (base$H6ub - base$L6ub),
                                      "M v. L by time interaction" = (base$M9ub - base$L9ub) - (base$M6ub - base$L6ub),
                                      "H v. M by fire interaction" = (base$H9b - base$M9b) - (base$H9ub - base$M9ub),
                                      "H v. L by fire interaction" = (base$H9b - base$L9b) - (base$H9ub - base$L9ub),
                                      "M v. L by fire interaction" = (base$M9b - base$L9b) - (base$M9ub - base$L9ub)), 
                               adjust = "mvt", type = "response")

cbind(summary(selected_contrasts),confint(selected_contrasts)[,5:6]) %>% 
  as.data.frame() %>% 
  mutate_at(2:8, round, 3)
