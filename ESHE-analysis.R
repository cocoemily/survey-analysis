#script for ESHE 
library(tidyverse)
library(MASS)
library(rstatix)
library(ggpubr)
library(ggthemes)
library(gapminder)
library(ggridges)
library(ggmosaic)
library(cowplot)
library(lmtest)
library(lme4)
library(vcd)
library(ggstatsplot)
library(rcompanion)
library(DescTools)
library(writexl)
library(dotwhisker)
library(OddsPlotty)
library(oddsratio)

source("site-comparison-functions.R")

theme_set(theme_bw())

artifacts1 = read_csv("data/cleaned_june_artifacts.csv")
artifacts2 = read_csv("data/cleaned_july_artifacts.csv")

paleocore_sss_artifact_form_all_versions_False_2022_08_01_05_36_18 <- read_delim("data/paleocore-sss_artifact-form_-_all_versions_-_False_-_2022-08-01-05-36-18.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

s10 = paleocore_sss_artifact_form_all_versions_False_2022_08_01_05_36_18 %>% filter(Site_name == "Semizbugu 10A")
s4 = paleocore_sss_artifact_form_all_versions_False_2022_08_01_05_36_18 %>% filter(Site_name == "Semizbugu 4")

artifacts = rbind(artifacts1, artifacts2)
collections = rbind(s10, s4)

artifacts = artifacts %>% mutate(location = ifelse(Site_name %in% c("Square 1", "Square 2", "Square 3"), "Semizbugu P1", 
                                                   ifelse(Site_name %in% c("Square 4", "Square 5"), "Semizbugu P2", "Semizbugu P5")))
artifacts$recycled = !is.na(artifacts$Recycling_description)
artifacts$double_patina = str_detect(artifacts$Recycling_indications, "double_patina")
artifacts$Raw_material_description = tolower(artifacts$Raw_material_description)

artifacts = artifacts %>% filter(is.na(Problem_notes))

collections$location = collections$Site_name
collections$recycled = !is.na(collections$Recycling_description)
collections$double_patina = str_detect(collections$Recycling_indications, "double_patina")
collections$Raw_material_description = tolower(collections$Raw_material_description)

collections = collections %>% filter(is.na(Problem_notes))

p1 = artifacts %>% filter(location == "Semizbugu P1")
p2 = artifacts %>% filter(location == "Semizbugu P2")
p5 = artifacts %>% filter(location == "Semizbugu P5")
s10a = collections %>% filter(location == 'Semizbugu 10A')
s4 = collections %>% filter(location == "Semizbugu 4")

all_artifacts = rbind(
  p1 %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Retouch, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight), 
  p2 %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Retouch, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight), 
  p5 %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Retouch, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight), 
  s10a %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Retouch, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight), 
  s4 %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Retouch, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight)
)

all_artifacts$Weathering_class = factor(all_artifacts$Weathering_class, 
                                        levels = c("strongly_weathered", "mildly_weathered", "weakly_weathered", "not_weathered", "other"))
all_artifacts$Artifact_type = factor(all_artifacts$Artifact_type, 
                                     levels = c("complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment", "shatter"))

all_artifacts = all_artifacts %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers


#LOCATION
ggplot(all_artifacts) +
  geom_bar(aes(recycled)) +
  facet_wrap(~location, ncol = 1) +
  coord_flip()

rl.table = table(all_artifacts %>% dplyr::select(location,recycled))
dim(rl.table)
dimnames(rl.table)
rl.table2 = apply(rl.table, c(1,2), sum)
rl.table3 = as.data.frame(rl.table2)
rl.table3$total = rl.table3$'FALSE' + rl.table3$'TRUE'
rl.table3$per_recycled = rl.table3$'TRUE'/rl.table3$total

chisq.test(rl.table2, simulate.p.value = T)
#different distributions between the locations

prop.test(x = rl.table2[,2], n = rl.table2[,2] + rl.table2[,1], alternative = "two.sided", correct = TRUE)
#at least one proportion of recycled objects at one site is significantly different than the other sites

pairwise.prop.test(x = rl.table2[,2], n = rl.table2[,2] + rl.table2[,1], alternative = "two.sided", correct = TRUE, p.adjust.method = "fdr")

test = cbind(rl.table2[,2], rl.table2[,1])
pairwise.prop.test(test, alternative = "two.sided", correct = TRUE, p.adjust.method = "bonferroni")


#### WEATHERING CLASS ####
ggplot(all_artifacts) +
  geom_bar(aes(Weathering_class, group = recycled, fill = recycled)) +
  facet_wrap(~location, ncol = 1) +
  coord_flip() +
  scale_fill_colorblind()

ggplot(all_artifacts) +
  geom_bar(aes(recycled)) +
  facet_wrap(~Weathering_class, ncol = 1) +
  coord_flip()
  

#weathering class distributions by location
wc.table = table(all_artifacts %>% dplyr::select(location, Weathering_class))
dim(wc.table)
wc.table2 = as.data.frame(ftable(wc.table))
wc.table2 = wc.table2 %>% group_by(location) %>%
  mutate(total = sum(Freq))
wc.table2 = wc.table2 %>%
  mutate(percent = Freq/total)

wide.wc = pivot_wider(wc.table2 %>% select(location, Weathering_class, percent), names_from = Weathering_class, values_from = percent)
write_xlsx(wide.wc, "~/Desktop//NYU/Dissertation-Research/papers/Coco_AK/figures/weathering.xlsx")

fisher.test(wc.table, simulate.p.value = T)
#distributions of weathering counts are not the same for each location

pairwiseNominalIndependence(wc.table, fisher = T, gtest= F, chisq = F, method = "fdr", digits = 3, simulate.p.value = T)
#no significant differences between P1 and P2? -- without bonferroni correction these two are significantly different
#no significant differences between P2 and P5


#### weathering class and recycling ####
wc.r.table = table(all_artifacts %>% dplyr::select(location, recycled, Weathering_class))
dim(wc.r.table)
dimnames(wc.r.table)
wc.r.table2 = apply(wc.r.table, c(1,2,3), sum)

wc.r = as.data.frame(ftable(wc.r.table))
Table = xtabs(Freq ~ location + recycled + Weathering_class,
              data=wc.r)
ftable(Table)

mantelhaen.test(Table)
groupwiseCMH(Table,
             group   = 1,
             fisher  = TRUE,
             gtest   = FALSE,
             chisq   = FALSE,
             method  = "fdr",
             correct = "none",
             digits  = 3, 
             simulate.p.value = T)
# Hypotheses
# •  Null hypothesis:  There is no association between the two inner variables.
# •  Alternative hypothesis (two-sided): There is an association between the two inner variables.
##suggests that there an association between recycling and weathering class at all sites

#regression analysis
reg.data = all_artifacts %>% select(location, recycled, Weathering_class)
logmodel = glm(recycled ~ Weathering_class, data = reg.data, family = binomial())
summary(logmodel)

logmodel2 = glmer(recycled ~ Weathering_class + (1 | location), data = reg.data, family = binomial)
print(logmodel2, corr = FALSE)
summary(logmodel2)


p1.glm = glm(recycled ~ Weathering_class, data = (reg.data %>% filter(location == "Semizbugu P1")), family = binomial)
summary(p1.glm)
p1.coef = tidy(p1.glm, conf.int = TRUE, exponentiate = T)
p1.coef$location = "Semizbugu P1"
p1.coef = p1.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
p2.glm = glm(recycled ~ Weathering_class, data = (reg.data %>% filter(location == "Semizbugu P2")), family = binomial)
summary(p2.glm)
p2.coef = tidy(p2.glm, conf.int = TRUE, exponentiate = T)
p2.coef$location = "Semizbugu P2"
p2.coef = p2.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
p5.glm = glm(recycled ~ Weathering_class, data = (reg.data %>% filter(location == "Semizbugu P5")), family = binomial)
summary(p5.glm)
p5.coef = tidy(p5.glm, conf.int = TRUE, exponentiate = T)
p5.coef$location = "Semizbugu P5"
p5.coef = p5.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
s10a.glm = glm(recycled ~ Weathering_class, data = (reg.data %>% filter(location == "Semizbugu 10A")), family = binomial)
summary(s10a.glm)
s10a.coef = tidy(s10a.glm, conf.int = TRUE, exponentiate = T)
s10a.coef$location = "Semizbugu 10A"
s10a.coef = s10a.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
s4.glm = glm(recycled ~ Weathering_class, data = (reg.data %>% filter(location == "Semizbugu 4")), family = binomial)
summary(s4.glm)
s4.coef = tidy(s4.glm, conf.int = TRUE, exponentiate = T)
s4.coef$location = "Semizbugu 4"
s4.coef = s4.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))


wc.coeffs = rbind(p1.coef, p5.coef, s4.coef, s10a.coef)
wc.coeffs$term = factor(wc.coeffs$term, levels = c("(Intercept)", "Weathering_classmildly_weathered", 
                                                   "Weathering_classweakly_weathered", 
                                                   "Weathering_classnot_weathered", 
                                                   "Weathering_classother"))
wc.coeffs = wc.coeffs %>%
  separate(location, c('site', "short_loc"), sep = " ")

wc.odds = ggplot(wc.coeffs, aes(x = term, y = estimate, shape = significance, color = significance))+
  geom_hline(yintercept = 1, color = "grey", linetype = 2) +
  geom_point(size = 3)+
  #geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  coord_flip() +
  facet_wrap(~short_loc, ncol = 1, strip.position = "right", 
             labeller = ) +
  scale_x_discrete(labels = c(
    "(Intercept)" = "(Intercept)", 
    "Weathering_classmildly_weathered" = "Mildly", 
    "Weathering_classweakly_weathered" = "Weakly", 
    "Weathering_classnot_weathered" = "Not", 
    "Weathering_classother" = "Other"
  )) +
  labs(x = expression(atop("Weathering category",atop("compared to strongly"))), y = "Odds ratio") +
  scale_color_colorblind() +
  theme(axis.title = element_text(size = 18), 
        legend.text = element_text(size = 14), legend.title = element_blank(), 
        axis.text = element_text(size = 14), strip.text = element_text(size = 14))
plot(wc.odds)

ggsave(filename = "~/Desktop/NYU/6th_Year/conferences/ESHE2022/poster-figures/weathering_recycling_odds.tiff", 
       plot = wc.odds, dpi = 300, width = 6, height = 4)


#### ARTIFACT TYPES ####
ggplot(all_artifacts) +
  geom_bar(aes(Artifact_type, group = recycled, fill = recycled)) +
  facet_wrap(~location, ncol = 1) +
  coord_flip() +
  scale_fill_colorblind()

atl.table = table(all_artifacts %>% dplyr::select(location,Artifact_type))
dim(atl.table)
dimnames(atl.table)
atl.table2 = apply(atl.table, c(1,2), sum)

#fishers exact test -- some small values
fisher.test(atl.table2, simulate.p.value = T)
#distributions of artifact types are not the same for each location

pairwiseNominalIndependence(atl.table2, fisher = T, gtest= F, chisq = F, method = "fdr", digits = 4, simulate.p.value = T)
#no significant differences between 10A and 4?
#no significant differences between 4 and P2?
#no significant differences between P1 and P2?
#no significant differences between P2 and P5

#### artifact type and recycling ####
at.r.table = table(all_artifacts %>% dplyr::select(location, recycled, Artifact_type))
dim(at.r.table)
dimnames(at.r.table)
at.r.table2 = apply(wc.r.table, c(1,2,3), sum)

at.r = as.data.frame(ftable(at.r.table))
Table = xtabs(Freq ~ location + recycled + Artifact_type,
              data=at.r)
ftable(Table)

mantelhaen.test(Table)
groupwiseCMH(Table,
             group   = 1,
             fisher  = TRUE,
             gtest   = FALSE,
             chisq   = FALSE,
             method  = "fdr",
             correct = "none",
             digits  = 3, 
             simulate.p.value = T)
# Hypotheses
# •  Null hypothesis:  There is no association between the two inner variables.
# •  Alternative hypothesis (two-sided): There is an association between the two inner variables.
##suggests that there an association between recycling and artifact type at all sites

#regression analysis
reg.data = all_artifacts %>% select(location, recycled, Artifact_type)
logmodel = glm(recycled ~ Artifact_type, data = reg.data, family = binomial())
summary(logmodel)

logmodel2 = glmer(recycled ~ Artifact_type + (1 | location), data = reg.data, family = binomial)
print(logmodel2, corr = FALSE)
summary(logmodel2)

p1.glm = glm(recycled ~ Artifact_type, data = (reg.data %>% filter(location == "Semizbugu P1")), family = binomial)
summary(p1.glm)
p1.coef = tidy(p1.glm, conf.int = TRUE, exponentiate = T)
p1.coef$location = "Semizbugu P1"
p1.coef = p1.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
p2.glm = glm(recycled ~ Artifact_type, data = (reg.data %>% filter(location == "Semizbugu P2")), family = binomial)
summary(p2.glm)
p2.coef = tidy(p2.glm, conf.int = TRUE, exponentiate = T)
p2.coef$location = "Semizbugu P2"
p2.coef = p2.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
p5.glm = glm(recycled ~ Artifact_type, data = (reg.data %>% filter(location == "Semizbugu P5")), family = binomial)
summary(p5.glm)
p5.coef = tidy(p5.glm, conf.int = TRUE, exponentiate = T)
p5.coef$location = "Semizbugu P5"
p5.coef = p5.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
s10a.glm = glm(recycled ~ Artifact_type, data = (reg.data %>% filter(location == "Semizbugu 10A")), family = binomial)
summary(s10a.glm)
s10a.coef = tidy(s10a.glm, conf.int = TRUE, exponentiate = T)
s10a.coef$location = "Semizbugu 10A"
s10a.coef = s10a.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
s4.glm = glm(recycled ~ Artifact_type, data = (reg.data %>% filter(location == "Semizbugu 4")), family = binomial)
summary(s4.glm)
s4.coef = tidy(s4.glm, conf.int = TRUE, exponentiate = T)
s4.coef$location = "Semizbugu 4"
s4.coef = s4.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))

at.coeffs = rbind(p1.coef, p2.coef, p5.coef, s4.coef, s10a.coef)
at.coeffs$term = factor(at.coeffs$term, levels = c(
  "(Intercept)", 
  "Artifact_typebroken_flake", 
  "Artifact_typetool", 
  "Artifact_typetool_fragment", 
  "Artifact_typecore", 
  "Artifact_typecore_fragment", 
  "Artifact_typeshatter"
))

at.coeffs = at.coeffs %>%
  separate(location, c('site', "short_loc"), sep = " ")

at.odds = ggplot(at.coeffs, aes(x = term, y = estimate, shape = significance, color = significance))+
  geom_hline(yintercept = 1, color = "grey", linetype = 2) +
  geom_point(size = 2.5)+
  #geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  coord_flip() +
  facet_wrap(~short_loc, ncol = 1, strip.position = "right") +
  scale_x_discrete(labels = c(
    "(Intercept)" = "(Intercept)", 
    "Artifact_typebroken_flake" = "Broken flake", 
    "Artifact_typetool" = "Tool", 
    "Artifact_typetool_fragment" = "Tool frag", 
    "Artifact_typecore" = "Core", 
    "Artifact_typecore_fragment" = "Core frag", 
    "Artifact_typeshatter" = "Shatter"
  )) +
  labs(x = expression(atop("Artifact type",atop("compared to complete flakes"))), y = "Odds ratio") +
  scale_color_colorblind() +
  theme(axis.title = element_text(size = 18), 
        legend.text = element_text(size = 14), legend.title = element_blank(), 
        axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 14))
plot(at.odds)

ggsave(filename = "~/Desktop/NYU/6th_Year/conferences/ESHE2022/poster-figures/artifact_type_recycling_odds.tiff", 
       plot = at.odds, dpi = 300, width = 6, height = 4)
  

#### retouch and recycling ####
r.r.table = table(all_artifacts %>% dplyr::select(location, recycled, Retouch))
dim(r.r.table)
dimnames(r.r.table)
r.r.table2 = apply(wc.r.table, c(1,2,3), sum)

r.r = as.data.frame(ftable(r.r.table))
Table = xtabs(Freq ~ location + recycled + Retouch,
              data=r.r)
ftable(Table)

mantelhaen.test(Table)
groupwiseCMH(Table,
             group   = 1,
             fisher  = TRUE,
             gtest   = FALSE,
             chisq   = FALSE,
             method  = "fdr",
             correct = "none",
             digits  = 3, 
             simulate.p.value = T)


#### RAW MATERIALS ####
ggplot(all_artifacts) +
  geom_bar(aes(Raw_material_description)) +
  facet_wrap(~location, ncol=1) +
  coord_flip()

unique(all_artifacts$Raw_material_description)

all_artifacts$grain_size = ifelse(
  str_detect(all_artifacts$Raw_material_description, "not"), "not fine grained", 
  ifelse(
    str_detect(all_artifacts$Raw_material_description, "less"), "less fine grained",
    "fine grained"
  )
)

unique(all_artifacts$grain_size)
ggplot(all_artifacts) +
  geom_bar(aes(grain_size)) +
  facet_wrap(~location, ncol = 1) +
  coord_flip()

gr.table = table(all_artifacts %>% dplyr::select(location, grain_size))
dim(gr.table)
dimnames(gr.table)
gr.table2 = apply(gr.table, c(1,2), sum)

pairwiseNominalIndependence(gr.table2, fisher = T, gtest= F, chisq = F, method = "fdr", digits = 4, simulate.p.value = T)

#### grain size and recycling ####
gr.r.table = table(all_artifacts %>% filter(!location %in% c("Semizbugu 10A", "Semizbugu 4")) %>% dplyr::select(location, recycled, grain_size))
dim(gr.r.table)
dimnames(gr.r.table)
gr.r.table2 = apply(wc.r.table, c(1,2,3), sum)

gr.r = as.data.frame(ftable(gr.r.table))
Table = xtabs(Freq ~ location + recycled + grain_size,
              data=gr.r)
ftable(Table)

mantelhaen.test(Table)
groupwiseCMH(Table,
             group   = 1,
             fisher  = TRUE,
             gtest   = FALSE,
             chisq   = FALSE,
             method  = "fdr",
             correct = "none",
             digits  = 3, 
             simulate.p.value = T)
# Hypotheses
# •  Null hypothesis:  There is no association between the two inner variables.
# •  Alternative hypothesis (two-sided): There is an association between the two inner variables.
##suggests that there an association between recycling and artifact type at all sites

#regression analysis
reg.data = all_artifacts %>% select(location, recycled, grain_size)
logmodel = glm(recycled ~ grain_size, data = reg.data, family = binomial())
summary(logmodel)


p1.glm = glm(recycled ~ grain_size, data = (reg.data %>% filter(location == "Semizbugu P1")), family = binomial)
summary(p1.glm)
p1.coef = tidy(p1.glm, conf.int = TRUE, exponentiate = T)
p1.coef$location = "Semizbugu P1"
p1.coef = p1.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
p2.glm = glm(recycled ~ grain_size, data = (reg.data %>% filter(location == "Semizbugu P2")), family = binomial)
summary(p2.glm)
p2.coef = tidy(p2.glm, conf.int = TRUE, exponentiate = T)
p2.coef$location = "Semizbugu P2"
p2.coef = p2.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))
p5.glm = glm(recycled ~ grain_size, data = (reg.data %>% filter(location == "Semizbugu P5")), family = binomial)
summary(p5.glm)
p5.coef = tidy(p5.glm, conf.int = TRUE, exponentiate = T)
p5.coef$location = "Semizbugu P5"
p5.coef = p5.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))


gs.coeffs = rbind(p1.coef, p2.coef, p5.coef)
gs.coeffs$term = factor(gs.coeffs$term, levels = c(
  "(Intercept)", 
  "grain_sizeless fine grained", 
  "grain_sizenot fine grained" 
  
))

gs.coeffs = gs.coeffs %>%
  separate(location, c("site", "short_loc"), sep = " ")

gs.odds = ggplot(gs.coeffs, aes(x = term, y = estimate, shape = significance, color = significance))+
  geom_hline(yintercept = 1, color = "grey", linetype = 2) +
  geom_point(size = 3)+
  #geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  coord_flip() +
  facet_wrap(~short_loc, ncol = 1, strip.position = "right") +
  scale_x_discrete(labels = c(
    "(Intercept)" = "(Intercept)",  
    "grain_sizeless fine grained" = "Less fine", 
    "grain_sizenot fine grained" = "Not fine" 
  )) +
  labs(x = expression(atop("Grain size",atop("compared to fine"))), y = "Odds ratio") +
  scale_color_colorblind() +
  theme(axis.title = element_text(size = 18), 
        legend.text = element_text(size = 14), legend.title = element_blank(), 
        axis.text = element_text(size = 14), strip.text = element_text(size = 14))
plot(gs.odds)

ggsave(filename = "~/Desktop/NYU/6th_Year/conferences/ESHE2022/poster-figures/grain_size_recycling_odds.tiff", 
       plot = gs.odds, dpi = 300, width = 6, height = 4)


### SIZE ####
ks.pair = all_artifacts %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 200 & Width <= 200 & Thickness <= 200) %>% #size of calipers
  select(Thickness, Length, Width, Weight, location)


## WEIGHT ####
#P1 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#P2 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#P5 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#10A comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#4 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Weight, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Weight, 
  alternative = "two.sided", 
  simulate.p.value = T
)


## LENGTH ####
#P1 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#P2 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#P5 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#10A comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#4 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Length, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Length, 
  alternative = "two.sided", 
  simulate.p.value = T
)


## WIDTH ####
#P1 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#P2 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#P5 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#10A comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#4 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Width, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Width, 
  alternative = "two.sided", 
  simulate.p.value = T
)


## THICKNESS ####
#P1 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P1"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#P2 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P2"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#P5 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu P5"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#10A comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 10A"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu 4"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

#4 comparisons
ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P2"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P5"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu 10A"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)

ks.test(
  x = (ks.pair %>% filter(location == "Semizbugu 4"))$Thickness, 
  y = (ks.pair %>% filter(location == "Semizbugu P1"))$Thickness, 
  alternative = "two.sided", 
  simulate.p.value = T
)


#### size and recycling ####
ks_test_recycled_vs_not(all_artifacts) 
reg.data = all_artifacts %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 200 & Width <= 200 & Thickness <= 200) %>%
  filter(Weight < 1000)

pd = reg.data %>% gather(key = "measurement", value = "value", Weight, Length, Width, Thickness) %>%
  separate(location, c("site", "short_loc", sep = " "))
pd$measurement = factor(pd$measurement, levels = c("Length", "Width", "Thickness", "Weight"))

stat.test1 = pd %>% group_by(measurement, short_loc) %>% t_test(value ~ recycled, ref.group = "FALSE") %>%
  add_xy_position(x = "recycled")

ks.df = data.frame(short_loc = c(
  rep("P1", 4),
  rep("P2", 4),
  rep("P5", 4),
  rep("10A", 4),
  rep("4", 4)
),
measurement = c(
  "Length", "Width", "Thickness", "Weight",
  "Length", "Width", "Thickness", "Weight",
  "Length", "Width", "Thickness", "Weight",
  "Length", "Width", "Thickness", "Weight",
  "Length", "Width", "Thickness", "Weight"
),
ks_p = c(
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P1"))$Length, (reg.data %>% filter(recycled != T & location == "Semizbugu P1"))$Length)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P1"))$Width, (reg.data %>% filter(recycled != T & location == "Semizbugu P1"))$Width)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P1"))$Thickness, (reg.data %>% filter(recycled != T & location == "Semizbugu P1"))$Thickness)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P1"))$Weight, (reg.data %>% filter(recycled != T & location == "Semizbugu P1"))$Weight)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P2"))$Length, (reg.data %>% filter(recycled != T & location == "Semizbugu P2"))$Length)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P2"))$Width, (reg.data %>% filter(recycled != T & location == "Semizbugu P2"))$Width)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P2"))$Thickness, (reg.data %>% filter(recycled != T & location == "Semizbugu P2"))$Thickness)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P2"))$Weight, (reg.data %>% filter(recycled != T & location == "Semizbugu P2"))$Weight)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P5"))$Length, (reg.data %>% filter(recycled != T & location == "Semizbugu P5"))$Length)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P5"))$Width, (reg.data %>% filter(recycled != T & location == "Semizbugu P5"))$Width)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P5"))$Thickness, (reg.data %>% filter(recycled != T & location == "Semizbugu P5"))$Thickness)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu P5"))$Weight, (reg.data %>% filter(recycled != T & location == "Semizbugu P5"))$Weight)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu 10A"))$Length, (reg.data %>% filter(recycled != T & location == "Semizbugu 10A"))$Length)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu 10A"))$Width, (reg.data %>% filter(recycled != T & location == "Semizbugu 10A"))$Width)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu 10A"))$Thickness, (reg.data %>% filter(recycled != T & location == "Semizbugu 10A"))$Thickness)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu 10A"))$Weight, (reg.data %>% filter(recycled != T & location == "Semizbugu 10A"))$Weight)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu 4"))$Length, (reg.data %>% filter(recycled != T & location == "Semizbugu 4"))$Length)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu 4"))$Width, (reg.data %>% filter(recycled != T & location == "Semizbugu 4"))$Width)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu 4"))$Thickness, (reg.data %>% filter(recycled != T & location == "Semizbugu 4"))$Thickness)$p),
  signif(ks.test((reg.data %>% filter(recycled == T & location == "Semizbugu 4"))$Weight, (reg.data %>% filter(recycled != T & location == "Semizbugu 4"))$Weight)$p)
))

stat.test2 = stat.test1 %>% left_join(ks.df, by = c("short_loc", "measurement"))

stat.test2 = stat.test2 %>% add_significance(
  p.col = "ks_p",
  output.col = "ks_p_sigf",
  cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
  symbols = c("****", "***", "**", "*", "ns")
)

pd$measurement = factor(pd$measurement, levels = c("Length", "Width", "Thickness", "Weight"))

plot = ggplot(pd, aes(x = as.factor(recycled), y = value, color = as.factor(recycled))) +
  geom_boxplot() +
  facet_grid(factor(measurement, levels = c("Length", "Width", "Thickness", "Weight")) ~ short_loc, scales = "free") +
  scale_x_discrete(labels = c(
    "FALSE" = "no", "TRUE" = "yes"
  )) +
  labs(x = "Recycled?", y = "Value", color = "Recycled?") +
  theme(legend.position = "none") +
  stat_pvalue_manual(stat.test2, label = "ks_p_sigf", tip.length = 0.01,
                     label.size = 5, bracket.nudge.y = -150, 
                     bracket.shorten = 0.5, bracket.size = 0.5, 
                     hide.ns = F) +
  scale_color_colorblind()  +
  theme(axis.title = element_text(size = 24), 
        legend.text = element_text(size = 24), legend.title = element_blank(), 
        axis.text = element_text(size = 18), strip.text = element_text(size = 22))

plot(plot)
ggsave(filename = "~/Desktop/NYU/6th_Year/conferences/ESHE2022/poster-figures/ks_recycled_vs_not_by_site.tiff", 
       plot = plot, dpi = 300, width = 8, height = 7)

all.log = glm(recycled ~ Length + Width + Thickness + Weight, data = reg.data, family = binomial)
summary(all.log)

dist_comp_recycled_vs_all(all_artifacts)

ks_test_recycled_vs_not(p1)
reg.data = p1 %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 175 & Width <= 175 & Thickness <= 175) #size of calipers
p1.log = glm(recycled ~ Length + Width + Thickness + Weight, data = reg.data, family = binomial)
summary(p1.log)

ks_test_recycled_vs_not(p2)
reg.data = p2 %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 175 & Width <= 175 & Thickness <= 175) #size of calipers
p2.log = glm(recycled ~ Length + Width + Thickness + Weight, data = reg.data, family = binomial)
summary(p2.log)

ks_test_recycled_vs_not(p5)
reg.data = p5 %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 175 & Width <= 175 & Thickness <= 175) #size of calipers
p5.log = glm(recycled ~ Length + Width + Thickness + Weight, data = reg.data, family = binomial)
summary(p5.log)

ks_test_recycled_vs_not(s10a)
reg.data = s10a %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 175 & Width <= 175 & Thickness <= 175) #size of calipers
s10a.log = glm(recycled ~ Length + Width + Thickness + Weight, data = reg.data, family = binomial)
summary(s10a.log)

ks_test_recycled_vs_not(s4)
reg.data = s4 %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 175 & Width <= 175 & Thickness <= 175) #size of calipers
s4.log = glm(recycled ~ Length + Width + Thickness + Weight, data = reg.data, family = binomial)
summary(s4.log)

recycled_dist_comp = function(data) {
  data = data %>%
    mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
           Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
           Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
    filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers
  
  pd = data %>% gather(key = "measurement", value = "value", Weight, Length, Width, Thickness)
  pd$measurement = factor(pd$measurement, levels = c("Length", "Width", "Thickness", "Weight"))
  
  p1 = ggplot(data = pd %>% filter(recycled == T)) + 
    geom_density(mapping = aes(x = log(value), 
                               color = recycled, fill = recycled), 
                 size = 0.5) + 
    geom_density(data = pd, mapping = aes(x = log(value)), 
                 color = "gray20", size = 0.5) + 
    facet_wrap(~measurement) +
    guides(color = "none", fill = "none") +
    scale_fill_manual(values = alpha("#56B4E9", 0.7)) + 
    scale_color_manual(values = alpha("#56B4E9", 1)) +
    ylim(0, 1) +
    xlim(-1, 8) +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=6,face="bold"), 
          title = element_text(size=8), 
          strip.text.x = element_text(size = 6))
  
  return(p1)
}

data = s4 %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers

summary(
  data %>% filter(recycled == T) %>% select(Length, Width, Thickness, Weight)
)

summary(
  data %>% filter(recycled == F) %>% select(Length, Width, Thickness, Weight)
)

#### REGRESSION ####
reg.data= all_artifacts %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers

logmodel = glm(recycled ~ Weathering_class + Artifact_type + grain_size + Thickness + Length + Width + Weight, data = reg.data, family = binomial())
summary(logmodel)

# logmodel2 = glmer(recycled ~ Weathering_class + Artifact_type + grain_size + Thickness + Length + Width + Weight + (1 | location), data = reg.data, family = binomial())
# print(logmodel2, corr = FALSE)
# summary(logmodel2)
# anova(logmodel2)

p1.rd = reg.data %>% filter(location == "Semizbugu P1")
p1.logmodel = glm(recycled ~ Weathering_class + Artifact_type + grain_size + Thickness + Length + Width + Weight, data = p1.rd, family = binomial())
summary(p1.logmodel)
p1.coef = tidy(p1.logmodel, conf.int = TRUE, exponentiate = T)
p1.coef$location = "Semizbugu P1"
p1.coef = p1.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))

p2.rd = reg.data %>% filter(location == "Semizbugu P2")
p2.logmodel = glm(recycled ~ Weathering_class + Artifact_type + grain_size + Thickness + Length + Width + Weight, data = p2.rd, family = binomial())
summary(p2.logmodel)
p2.coef = tidy(p2.logmodel, conf.int = TRUE, exponentiate = T)
p2.coef$location = "Semizbugu P2"
p2.coef = p2.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))

p5.rd = reg.data %>% filter(location == "Semizbugu P5")
p5.logmodel = glm(recycled ~ Weathering_class + Artifact_type + grain_size + Thickness + Length + Width + Weight, data = p5.rd, family = binomial())
summary(p5.logmodel)
p5.coef = tidy(p5.logmodel, conf.int = TRUE, exponentiate = T)
p5.coef$location = "Semizbugu P5"
p5.coef = p5.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))

s4.rd = reg.data %>% filter(location == "Semizbugu 4") %>%
  mutate(recycl = ifelse(recycled, 1, 0))
s4.logmodel = glm(recycled ~ Weathering_class + Artifact_type + grain_size + Thickness + Length + Width + Weight, data = s4.rd, family = binomial())
summary(s4.logmodel)
s4.coef = tidy(s4.logmodel, conf.int = TRUE, exponentiate = T)
s4.coef$location = "Semizbugu 4"
s4.coef = s4.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))

s10a.rd = reg.data %>% filter(location == "Semizbugu 10A") %>%
  mutate(recycl = ifelse(recycled, 1, 0))
s10a.logmodel = glm(recycled ~ Weathering_class + Artifact_type  + Thickness + Length + Width + Weight, data = s10a.rd, family = binomial())
summary(s10a.logmodel)

s10a.coef = tidy(s10a.logmodel, conf.int = TRUE, exponentiate = T)
s10a.coef$location = "Semizbugu 10A"
s10a.coef = s10a.coef %>% 
  mutate(significance = ifelse(p.value <= 0.05, "signif.", "not signif."))


# s10a.logmodel2 = caret::train(
#   recycl ~ Weathering_class + Artifact_type  + Thickness + Length + Width + Weight,
#   data = s10a.rd,
#   method = "glm",
#   family = "binomial",
#   na.action = na.omit
# )
# summary(s10a.logmodel2)
# odds_plot(s10a.logmodel2$finalModel)


coeffs = rbind(p1.coef, p2.coef, p5.coef, s4.coef, s10a.coef)

ggplot(coeffs, aes(x = term, y = estimate, shape = significance, color = significance))+
  geom_hline(yintercept = 1, color = "grey", linetype = 2) +
  geom_point()+
  #geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  coord_flip() +
  facet_grid(~location, scales = "free" )

gs.coeff = coeffs %>% filter(str_detect(term, "grain_size"))
at.coeff = coeffs %>% filter(str_detect(term, "Artifact_type"))

ggplot(at.coeff, aes(x = term, y = estimate, shape = significance, color = significance))+
  geom_hline(yintercept = 1, color = "grey", linetype = 2) +
  geom_point()+
  #geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  coord_flip() +
  facet_grid(~location, scales = "free" )

#### MOSAIC PLOTS ####
all_artifacts$Weathering_class = factor(all_artifacts$Weathering_class, 
                                        levels = c("strongly_weathered", "mildly_weathered", "weakly_weathered", "not_weathered", "other"))
wm = ggplot(data = all_artifacts) +
  geom_mosaic(aes(x = product(Weathering_class, location), fill = Weathering_class), na.rm = T, offset = 0.04) +
  scale_fill_colorblind(labels = c(
    "Strongly weathered", "Mildly weathered", "Weakly weathered", "Not weathered", "Other"
    )) +
  scale_y_productlist(name = "Weathering category",
    labels = c(
    "Strongly", "Mildly ", "Weakly", "Not", "Other"
  )) +
  scale_x_productlist(name = "Location",
    labels = c(
    "4", "10A", "P1", "P2", "P5"
  )) +
  #guides(fill = "Weathering category") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14), strip.text = element_text(size = 12), 
        legend.position = "none") 
plot(wm)
ggsave(filename = "~/Desktop/NYU/6th_Year/conferences/ESHE2022/poster-figures/weathering_mosaic.tiff", 
       plot = wm, dpi = 300, width = 4, height = 4)
  

atm = ggplot(data = all_artifacts) +
  geom_mosaic(aes(x = product(Artifact_type, location), fill = Artifact_type), na.rm = T, offset = 0.04) +
  scale_fill_colorblind(labels = c(
    "Complete flake", "Broken flake", "Tool", "Tool frag", "Core", "Core frag", "Shatter"
  )) +
  scale_y_productlist(name = "Artifact type",
                      labels = c(
                        "Complete flake", "Broken flake", "Tool", "Tool frag", "Core", "Core frag", "Shatter"
                      )) +
  scale_x_productlist(name = "Location",
                      labels = c(
                        "4", "10A", "P1", "P2", "P5"
                      )) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14), strip.text = element_text(size = 12), 
        legend.position = "none")
plot(atm)
ggsave(filename = "~/Desktop/NYU/6th_Year/conferences/ESHE2022/poster-figures/artifact_type_mosaic.tiff", 
       plot = atm, dpi = 300, width = 4, height = 4)

gsm = ggplot(data = all_artifacts) +
  geom_mosaic(aes(x = product(grain_size, location), fill = grain_size), na.rm = T, offset = 0.04) +
  scale_fill_colorblind(labels = c(
    "Fine", "Less fine", "Not fine"
  )) +
  scale_y_productlist(name = "Grain size",
                      labels = c(
                        "Fine", "Less fine", "Not fine"
                      )) +
  scale_x_productlist(name = "Location",
                      labels = c(
                        "4", "10A", "P1", "P2", "P5"
                      )) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14), strip.text = element_text(size = 12), 
        legend.position = "none")
ggsave(filename = "~/Desktop/NYU/6th_Year/conferences/ESHE2022/poster-figures/grain_size_mosaic.tiff", 
       plot = gsm, dpi = 300, width = 4, height = 4)


#### DISTRIBUTION PLOTS ####
pd = ks.pair %>% gather(key = "measurement", value = "value", Weight, Length, Width, Thickness) %>%
  separate(location, c("site", "short_loc", sep = " "))
pd$measurement = factor(pd$measurement, levels = c("Length", "Width", "Thickness", "Weight"))


dist = ggplot(data = pd) + 
  # geom_density(data = pd, mapping = aes(x = log(value)), 
  #              color = "gray20", size = 0.5) + 
  geom_density(mapping = aes(x = log(value), 
                             fill = measurement, color = measurement), 
               size = 0.5) + 
  facet_grid(factor(measurement, levels = c("Length", "Width", "Thickness", "Weight")) ~short_loc) +
  guides(color = "none", fill = "none") +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  scale_x_continuous(breaks = c(0, 2, 5, 8)) +
  theme(axis.title = element_text(size = 24), 
        legend.text = element_text(size = 24), legend.title = element_blank(), 
        axis.text = element_text(size = 18), strip.text = element_text(size = 22))
plot(dist)

ggsave(filename = "~/Desktop/NYU/6th_Year/conferences/ESHE2022/poster-figures/size_distributions.tiff", 
       plot = dist, dpi = 300, width = 8, height = 7)
