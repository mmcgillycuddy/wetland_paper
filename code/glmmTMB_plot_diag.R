library(glmmTMB)
library(tidyverse)
library(ggplot2)
library(sjPlot)

# read in fire_water_model
load("results/comp_water_fire_March22.RData")
### --------------------------------------
### Plot random effects - diag
### --------------------------------------
object <- comp_water_fire_fit
listname <- "cond"
cnms <- object$modelInfo$reTrms[[listname]]$cnms   ## list of (named) terms and X columns
reStruc <- object$modelInfo$reStruc[[paste0(listname, "ReStruc")]] ## random-effects structure
flist <- object$modelInfo$reTrms[[listname]]$flist ## list of grouping variables
levs <- lapply(flist, levels)

pl <- object$obj$env$parList(object$fit$par, object$fit$parfull)
nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
nb <- vapply(reStruc, function(x) x$blockReps, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)
bc <- vapply(reStruc, function(x) x$blockCode, numeric(1)) ## code block
nbseq <- rep.int(seq_along(nb), nb * nc)       ## splitting vector
ml.b <- split(pl$b, nbseq)
ml <- ml.b
for (i in seq_along(ml.b)) {
  ml[[i]] <- matrix(ml.b[[i]], ncol = nc[i], byrow = TRUE,
                    dimnames = list(NULL, cnms[[i]]))
}

block.plot <- 1 # diag re is the first random effect in model
b <- as.matrix(ml[[block.plot]])
rownames(b) <- levs[[block.plot]]
b.plot <- t(b)
b.plot.long <- as.data.frame(b.plot) %>%
  tibble::rownames_to_column("Coefficient") %>%
  tidyr::pivot_longer(cols = -Coefficient,
                      names_to = c(names(levs[block.plot])),
                      values_to = "b" )

s1 <- TMB::sdreport(object$obj,  getJointPrecision = TRUE)
# subset first
s1.joint <- s1$jointPrecision
s1.b <- s1.joint[rownames(s1.joint)=="b", colnames(s1.joint)=="b"]
b.cols <- split(1:ncol(s1.b), nbseq)
b.cols.plot <-  b.cols[[block.plot]]
s1.b.plot <- s1.b[b.cols.plot, b.cols.plot]
s1.b.inv <- solve(s1.b.plot)

diag.H.b <- as.data.frame(diag(s1.b.inv))
coef.nam <- rep(cnms[[block.plot]], nb[[block.plot]])
nam.grp <- rep(levs[[block.plot]], each =nc[[block.plot]])
var.re <- cbind(coef.nam, nam.grp, diag.H.b)

names(var.re) <- c("Coefficient", names(levs[block.plot]), "var.u")
save.image("comp_fire_water_plot_dataMarch22.RData")
load("results/comp_water_fire_plot_dataMarch22.RData")
ranef.plot <- left_join(b.plot.long, var.re, by = intersect(names(b.plot.long), names(var.re)))
ranef.plot$sd.u <- sqrt(ranef.plot$var.u)

highlight_coef = c("Time_pt:water_treatmentL",
                   # "Time_pt:water_treatmentM",
                   "water_treatmentL:fire_treatmentb", 
                   "water_treatmentM:fire_treatmentb")
# add fixed effects
fixed_coefs <- fixef(comp_water_fire_fit)$cond %>% 
  as.data.frame() %>%
  rownames_to_column( var = "Coefficient") %>% 
  rename(fixed_coef = 2)

ranef_fixed_plot <- ranef.plot %>%
  mutate(conf.low = b - 1.96*sd.u,
         conf.high =  b + 1.96*sd.u,
         not.zero = ifelse(conf.high < 0 | conf.low > 0, 1, 0),
         comp_high = ifelse(Coefficient %in% highlight_coef &  not.zero, 1, 0 ) ) %>% 
  left_join(fixed_coefs) %>% 
  mutate(btot = b + fixed_coef,
         conf.low_tot = conf.low + fixed_coef,
         conf.high_tot = conf.high + fixed_coef,
         not.zero_tot = ifelse(conf.high_tot < 0 | conf.low_tot > 0, 1, 0),
         abund_high = ifelse(Coefficient %in% highlight_coef &  not.zero_tot, 1, 0 ),
         effect = ifelse(abund_high==1 & comp_high == 1,"abundance and composition",
                         ifelse(abund_high==1 & comp_high == 0,"abundance only",
                                ifelse(comp_high == 1, "composition only","neither"))))
# abundance = factor(highlight_tot, levels = c("1","0"), labels  = c("effect","no effect" )),
# composition = factor(highlight, levels = c("1","0"), labels  = c("effect", "no effect" )))
cols <- c("1" = "red", "0" = "black")
plot.spp <- ranef_fixed_plot$Species[ranef_fixed_plot$effect!="neither"] 

# New facet label names for supp variable
supp.labs <- c("L-H x time","(b-ub) x (L-H)","(b-ub) x (M-H)")
names(supp.labs) <- highlight_coef

ranef_fixed_plot %>% 
  filter(Species %in% plot.spp,
         Coefficient %in% highlight_coef,
         effect!="neither") %>% 
  ggplot(aes(x = as.factor(Species),
             y = btot,
             color = effect,
             linetype = effect)) +
  geom_point()+ 
  geom_linerange(aes(ymin = conf.low_tot, ymax = conf.high_tot)) +
  geom_hline(yintercept = 0) +
  geom_hline(data = fixed_coefs %>% filter(Coefficient %in% highlight_coef), 
             aes(yintercept = fixed_coef), linetype  = "dashed", color = "gray60")+ 
  facet_grid(~ Coefficient,
             scales = "free",
             labeller = labeller(Coefficient = supp.labs),
             switch = "y") +
  coord_flip() + theme_classic() +
  xlab("Species") + ylab("OR") +
  theme( axis.text = element_text( size = 12 ),
         axis.text.y = element_text( size = 10),
         axis.text.x = element_text( size = 8),
         axis.title = element_text( size = 12),
         strip.text = element_text(size = 8),
         legend.position="bottom",
         strip.text.x = element_text(size=8, angle=0)) +
  scale_color_manual(values  = c("black","black","gray60"))+
  scale_linetype_manual(values = c("dashed","solid","dashed"))

### PDirection

# Time_pt:water_treatmentL - 

# save.image("results/comp_fire_water_plot_dataMarch22.RData")

nd = PA_data %>% 
  group_by(Species,Swamp, veg, Time_pt,water_treatment,fire_treatment) %>% 
  count() %>% 
  mutate(Sod = "NA") %>% 
  filter(Species %in% plot.spp,
         veg == "tt",
         Swamp == "CA")

nd$pred = predict(comp_water_fire_fit, newdata = nd, se.fit = F, type = "response", re.form = NULL, allow.new.levels = TRUE)

save(nd,file = "pred_no_se.RData")

for(sp in plot.spp){
  p = nd %>%
    filter(Species == sp) %>% 
    ggplot(aes(Time_pt,pred, color = water_treatment, linetype = fire_treatment))+
    geom_point()+
    geom_line()+
    theme_classic()
  ggsave(filename = paste0("results/spplots/plot_",sp,".png"), device = "png", plot  = p)
}

