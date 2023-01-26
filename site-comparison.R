#comparison of different sites
library(tidyverse)
library(MASS)
library(rstatix)
library(ggpubr)
library(ggthemes)
library(gapminder)
library(ggridges)
library(cowplot)
library(lmtest)
library(vcd)
library(ggstatsplot)
library(rcompanion)
library(lme4)
library(QuantPsyc)

source("site-comparison-functions.R")

theme_set(theme_bw())

#### DATA LOADING ####
artifacts1 = read_csv("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/cleaned_june_artifacts.csv")
artifacts2 = read_csv("~/Desktop/NYU/Dissertation-Research/Survey/July-survey/cleaned_july_artifacts.csv")

paleocore_sss_artifact_form_all_versions_False_2022_08_01_05_36_18 <- read_delim("~/Desktop/NYU/Dissertation-Research/Survey/July-survey/paleocore-sss_artifact-form_-_all_versions_-_False_-_2022-08-01-05-36-18.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

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

cols = c("Id_number", "location", "recycled", "double_patina", "Raw_material_description", 
         "Weathering_class", "Artifact_type", "Bordian_type", "Tool_type", "Flake_type", "Blank_form",
         "Dorsal_flake_scar_count", "Cortex_percentage", "Flake_fragment", "Flake_termination",
         "Retouch", "Retouch_side", "Edge_damage",
         "Platform_thickness", "Platform_width",
         "Flake_thickness", "Flake_length", "Flake_width", 
         "Maximum_core_length", "Maximum_core_width", "Maximum_core_thickness", 
         "Weight")

all_artifacts = rbind(
  p1 %>% select_at(cols), 
  p2 %>% select_at(cols), 
  p5 %>% select_at(cols), 
  s10a %>% select_at(cols), 
  s4 %>% select_at(cols)
  )


#### DATA CLEANING ####
all_artifacts$Weathering_class = factor(all_artifacts$Weathering_class, 
                                        levels = c("strongly_weathered", "mildly_weathered", "weakly_weathered", "not_weathered", "other"))
all_artifacts$Artifact_type = factor(all_artifacts$Artifact_type, 
                                     levels = c("complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment", "shatter"))
all_artifacts$Flake_termination = factor(all_artifacts$Flake_termination, 
                                         levels = c("feather", "hinge", "plunge", "step", "other"))
all_artifacts$Flake_fragment = factor(all_artifacts$Flake_fragment , 
                                         levels = c("proximal", "medial", "distal", "other"))

all_artifacts$retouch.side = ifelse(str_detect(all_artifacts$Retouch_side, pattern = " "), "bifacial", all_artifacts$Retouch_side)
all_artifacts$retouch.side = factor(all_artifacts$retouch.side, levels = c("dorsal", "ventral", "bifacial"))
#all_artifacts$Retouch_side = factor(all_artifacts$Retouch_side, levels = c("dorsal", "ventral", "bifacial"))

all_artifacts$tool.type = ifelse(str_detect(all_artifacts$Tool_type, "notch denticulate"), "notch/denticulate", 
                        ifelse(str_detect(all_artifacts$Tool_type, " "), "multiple", 
                               all_artifacts$Tool_type))
all_artifacts$tool.type = factor(all_artifacts$tool.type, levels = c("notch", "denticulate", "notch/denticulate", "scraper", "point", "biface", "multiple", "other"))

all_artifacts = all_artifacts %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers

all_artifacts = all_artifacts %>%
  mutate(flake.type = ifelse(is.na(Blank_form), Flake_type, Blank_form))

source("bordian_types_dictionary.R")
all_artifacts$Bordian_name = ""
for(i in 1:nrow(all_artifacts)) {
  if(!is.na(all_artifacts$Bordian_type[i])) {
    all_artifacts$Bordian_name[i] = bordian_types[all_artifacts$Bordian_type[i]]
  }
}
all_artifacts$Bordian_name = factor(all_artifacts$Bordian_name, levels = bordian_levels)

all_artifacts = all_artifacts %>%
  mutate(
    Tool_type = ifelse(str_detect(Tool_type, "notch denticulate"), "notch/denticulate", 
                       ifelse(str_detect(Tool_type, " "), "multiple", 
                              Tool_type)), 
    Flake_type = ifelse(str_detect(Flake_type, pattern = "flake"), "flake", Flake_type)
    )
all_artifacts$Tool_type = factor(all_artifacts$Tool_type, levels = c("notch", "denticulate", "notch/denticulate", "scraper", "point", "biface", "multiple", "other"))
all_artifacts$Flake_type = factor(all_artifacts$Flake_type, levels = c("flake", "blade", "bladelet", "other"))


#### TO DO: raw material cleaning ####

#possible artifact rolling
all_artifacts$poss_roll = ifelse(all_artifacts$Bordian_name == "alternate scraper", TRUE, FALSE)
all_artifacts$poss_roll = ifelse(str_detect(all_artifacts$Retouch_side, pattern = " "), TRUE, all_artifacts$poss_roll)

rolled.rcycl = all_artifacts %>% filter(poss_roll == T & recycled == T)
all_artifacts = subset(all_artifacts, !(Id_number %in% rolled.rcycl$Id_number))


#### cleaning up environment ####
rm(list = c("artifacts", "artifacts1", "artifacts2", "collections", "paleocore_sss_artifact_form_all_versions_False_2022_08_01_05_36_18", "s10", "rolled.rcycl"))


#### Differences by location ####
ggplot(all_artifacts) +
  geom_bar(aes(recycled)) +
  facet_wrap(~location) +
  coord_flip()

rl.table = table(all_artifacts %>% dplyr::select(location,recycled))
chisq.test(rl.table)
pairwiseNominalIndependence(rl.table, simulate.p.value = T, 
                            fisher = T, chisq = F, gtest = F)
#No difference in distribution of recycled vs non recycled objects between P2 and P5

at.table = table(all_artifacts %>% dplyr::select(location, Artifact_type))
pairwiseNominalIndependence(at.table, simulate.p.value = T, 
                            fisher = T, chisq = F, gtest = F)
#no difference between S10A and S4 or between P2 and P5 for artifact type distributions

at.r.table = table(all_artifacts %>% dplyr::select(location, recycled, Artifact_type))
at.r.table2 = apply(at.r.table, c(1,2,3), sum)
at.r = as.data.frame(ftable(at.r.table))
Table = xtabs(Freq ~ location + recycled + Artifact_type, data=at.r)
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
##suggests that there an association between recycling and artifact type at all sites, except for Semizbugu 10A

ggplot(all_artifacts) +
  geom_bar(aes(x = Artifact_type, fill = recycled)) +
  facet_wrap(~location) +
  coord_flip() +
  scale_fill_colorblind()

##error where model is nearly unidentifiable -- likely due to low numbers of certain artifact types
rafit = glmer(recycled ~ Artifact_type + (1 | location), family = binomial(), data = all_artifacts)
summary(rafit)



wc.table = table(all_artifacts %>% dplyr::select(location, Weathering_class))
pairwiseNominalIndependence(wc.table, simulate.p.value = T, 
                            fisher = T, chisq = F, gtest = F)
#no difference between P2 and P5 in distribution of artifacts among weathering classes

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, Retouch)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)
#Semizbugu 4 has significant differences in distribution of retouch (presence/absence) compared to all other sites

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, Edge_damage)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, retouch.side)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)

rs.r.table = table(all_artifacts %>% dplyr::select(location, recycled, retouch.side))
rs.r.table2 = apply(rs.r.table, c(1,2,3), sum)
rs.r = as.data.frame(ftable(rs.r.table))
Table = xtabs(Freq ~ location + recycled + retouch.side,
              data=rs.r)
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
##no association between recycling and retouch side at any sites


pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, Bordian_name)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, tool.type)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)
#no difference in tool type distributions between 10A and 4, between P1 and P2, or between P2 and P5

tt.r.table = table(all_artifacts %>% dplyr::select(location, recycled, tool.type))
tt.r.table2 = apply(tt.r.table, c(1,2,3), sum)
tt.r = as.data.frame(ftable(tt.r.table))
Table = xtabs(Freq ~ location + recycled + tool.type,
              data=tt.r)
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
##suggests that there an association between recycling and tool type at P5 and S4, but not at the other sites

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, flake.type)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, Flake_fragment)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, Flake_termination)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)


pairwiseKS(all_artifacts %>% dplyr::select(location, Weight))

pairwiseKS(all_artifacts %>% dplyr::select(location, Length))

pairwiseKS(all_artifacts %>% dplyr::select(location, Width))

pairwiseKS(all_artifacts %>% dplyr::select(location, Thickness))

pairwiseKS(all_artifacts %>% dplyr::select(location, Cortex_percentage))

pairwiseKS(all_artifacts %>% dplyr::select(location, Platform_thickness))
pairwiseKS(all_artifacts %>% dplyr::select(location, Platform_width))

pairwiseKS(all_artifacts %>% dplyr::select(location, Dorsal_flake_scar_count))


#### Regressions between locations ####
##does location explain more variation than other variables?
reg.data = all_artifacts %>% select("recycled", "Weathering_class", 
                                    "Dorsal_flake_scar_count", "Cortex_percentage",
                                    "Thickness", "Length", "Width", 
                                    "Weight")
reg.data$recycled = as.factor(reg.data$recycled)

fit1 = glm(recycled ~ ., family = binomial(), data = reg.data)
summary(fit1)

reg.data2 = all_artifacts %>% select("location", "recycled", "Weathering_class", 
                                     "Dorsal_flake_scar_count", "Cortex_percentage",
                                     "Thickness", "Length", "Width", 
                                     "Weight")
fit2 = glm(recycled ~ ., family = binomial(), data = reg.data2)
summary(fit2)

anova(fit1, fit2, test = "LR")
##need to include location, but weakly and not weathered artifacts are still less recycled, 
##and overall long and wide artifacts are more likely to be recycled, but also thinner artifacts

###### standardized coefficients#####
scoef = as.data.frame(lm.beta::lm.beta(fit2)[["standardized.coefficients"]])
scoef$var = rownames(scoef)
colnames(scoef) = c("coef", "var")
rownames(scoef) = NULL
scoef$abs_coef = abs(scoef$coef)
scoef = scoef %>% filter(!is.na(abs_coef))

scoef$var = factor(scoef$var, 
                   levels = c(
                     "locationSemizbugu 4", "locationSemizbugu P1", "locationSemizbugu P2", "locationSemizbugu P5", 
                     "Weathering_classmildly_weathered", "Weathering_classweakly_weathered", "Weathering_classnot_weathered", "Weathering_classother",
                     "Dorsal_flake_scar_count", "Cortex_percentage", 
                     "Length", "Width", "Thickness", "Weight"
                   ))
scoef$var = fct_rev(scoef$var)

ggplot(scoef) +
  geom_hline(yintercept = 0, color = I("red")) +
  geom_hline(yintercept = min((scoef %>% filter(str_detect(var, "location")))$abs_coef), color = I("grey"), linetype = "dashed") +
  geom_point(aes(x = var, y = abs_coef, color = var))+
  coord_flip() +
  theme(legend.position = "none")


#####multiple regressions for each of the locations #####
reg.data.p1 = all_artifacts %>% filter(location == "Semizbugu P1") %>%
  select("recycled", "Weathering_class",
         "Dorsal_flake_scar_count", "Cortex_percentage",
         "Length", "Width", "Thickness", 
         "Weight")

fit.p1 = glm(recycled ~ ., family = binomial(), data = reg.data.p1)
summary(fit.p1) 
p1.df = tidy(fit.p1)
p1.df$location = "Semizbugu P1"
p1.df$signf = ifelse(p1.df$p.value < 0.05, "signf", "not signf")



reg.data.p2 = all_artifacts %>% filter(location == "Semizbugu P2") %>%
  select("recycled", "Weathering_class",
         "Dorsal_flake_scar_count", "Cortex_percentage",
         "Length", "Width", "Thickness", 
         "Weight")

fit.p2 = glm(recycled ~ ., family = binomial(), data = reg.data.p2)
summary(fit.p2) 
p2.df = tidy(fit.p2)
p2.df$location = "Semizbugu P2"
p2.df$signf = ifelse(p2.df$p.value < 0.05, "signf", "not signf")

reg.data.p5 = all_artifacts %>% filter(location == "Semizbugu P5") %>%
  select("recycled", "Weathering_class",
         "Dorsal_flake_scar_count", "Cortex_percentage",
         "Length", "Width", "Thickness", 
         "Weight")

fit.p5 = glm(recycled ~ ., family = binomial(), data = reg.data.p5)
summary(fit.p5) 
p5.df = tidy(fit.p1)
p5.df$location = "Semizbugu P5"
p5.df$signf = ifelse(p5.df$p.value < 0.05, "signf", "not signf")

reg.data.s10a = all_artifacts %>% filter(location == "Semizbugu 10A") %>%
  select("recycled", "Weathering_class",
         "Dorsal_flake_scar_count", "Cortex_percentage",
         "Length", "Width", "Thickness", 
         "Weight")

fit.s10a = glm(recycled ~ ., family = binomial(), data = reg.data.s10a)
summary(fit.s10a) 
s10a.df = tidy(fit.s10a)
s10a.df$location = "Semizbugu 10A"
s10a.df$signf = ifelse(s10a.df$p.value < 0.05, "signf", "not signf")

reg.data.s4 = all_artifacts %>% filter(location == "Semizbugu 4") %>%
  select("recycled", "Weathering_class",
         "Dorsal_flake_scar_count", "Cortex_percentage",
         "Length", "Width", "Thickness", 
         "Weight")

fit.s4 = glm(recycled ~ ., family = binomial(), data = reg.data.s4)
summary(fit.s4) 
s4.df = tidy(fit.s4)
s4.df$location = "Semizbugu 4"
s4.df$signf = ifelse(s4.df$p.value < 0.05, "signf", "not signf")

#jtools::plot_summs(fit.p1, fit.p2, fit.p5, fit.s10a, fit.s4)

allregs = rbind(p1.df, p2.df, p5.df, s10a.df, s4.df)
allregs$lower = allregs$estimate - allregs$std.error
allregs$upper = allregs$estimate + allregs$std.error

ggplot(allregs %>% filter(signf == "signf")) +
  geom_hline(yintercept = 0, color = I("red")) +
  geom_linerange(aes(x=term, ymin = lower, ymax = upper, group = location, color = location), position = position_dodge2(width = 0.5)) +
  geom_point(aes(x=term, y = estimate, group = location, color = location, shape = location), position = position_dodge2(width = 0.5), size = 2) +
  coord_flip() +
  scale_color_colorblind() +
  scale_shape_manual(values = c(0:2,5))

