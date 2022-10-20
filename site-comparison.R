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

source("site-comparison-functions.R")

theme_set(theme_bw())

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
         "Weathering_class", "Artifact_type", "Bordian_type", "Tool_type", "Flake_type", 
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

all_artifacts$Weathering_class = factor(all_artifacts$Weathering_class, 
                                        levels = c("strongly_weathered", "mildly_weathered", "weakly_weathered", "not_weathered", "other"))
all_artifacts$Artifact_type = factor(all_artifacts$Artifact_type, 
                                     levels = c("complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment", "shatter"))
all_artifacts$Flake_termination = factor(all_artifacts$Flake_termination, 
                                         levels = c("feather", "hinge", "plunge", "step", "other"))
all_artifacts$Flake_fragment = factor(all_artifacts$Flake_fragment , 
                                         levels = c("proximal", "medial", "distal", "other"))

all_artifacts = all_artifacts %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers

all_artifacts = all_artifacts %>%
  mutate(Retouch_side = ifelse(str_detect(all_artifacts$Retouch_side, pattern = " "), "bifacial", Retouch_side)) %>%
  mutate()




#### Differences by location ####
ggplot(all_artifacts) +
  geom_bar(aes(recycled)) +
  facet_wrap(~location) +
  coord_flip()

rl.table = table(all_artifacts %>% dplyr::select(location,recycled))
chisq.test(rl.table)
pairwiseNominalIndependence(rl.table)
#No difference in distribution of recycled vs non recycled objects between P2 and P5

at.table = table(all_artifacts %>% dplyr::select(location, Artifact_type))
at.table
pairwiseNominalIndependence(at.table, simulate.p.value = T, 
                            fisher = T, chisq = F, gtest = F)

wc.table = table(all_artifacts %>% dplyr::select(location, Weathering_class))
wc.table
pairwiseNominalIndependence(wc.table, simulate.p.value = T, 
                            fisher = T, chisq = F, gtest = F)

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, Retouch)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, Edge_damage)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)

pairwiseNominalIndependence(
  table(all_artifacts %>% dplyr::select(location, Retouch_side)),
  simulate.p.value = T, 
  fisher = T, chisq = F, gtest = F
)

#### old code ####
rl.table = table(artifacts %>% dplyr::select(location,recycled))
dim(rl.table)
dimnames(rl.table)
rl.table2 = apply(rl.table, c(1,2), sum)

M3 = loglm( ~ recycled * location, dat=rl.table2, fit=TRUE)
M4 = update(M3, ~ . - location:recycled)
#anova(M3, M4)
chisq.test(rl.table2)
#results indicate that there is a significant interaction between location and recycling
#this means that there is a difference in the recycled object amounts between locations


#Artifact types by location
ggplot(artifacts) +
  geom_bar(aes(Artifact_type)) +
  facet_wrap(~location) +
  coord_flip()

atl.table = table(artifacts %>% dplyr::select(location,Artifact_type))
dim(atl.table)
dimnames(atl.table)
atl.table2 = apply(atl.table, c(1,2), sum)

M5 = loglm( ~ Artifact_type * location, dat=atl.table2, fit=TRUE)
M6 = update(M5, ~ . - location:Artifact_type)
#anova(M5, M6)
chisq.test(atl.table2)
#results indicate that there is a significant interaction between location and artifact type
#this means that there is a difference in the types of artifact amounts between locations

  

## what types of objects are recycled?
recycled_objects_artifact_type = function(data) {
  df1 = data %>% group_by(location, Artifact_type, recycled) %>%
    summarize(count = n()) %>%
    group_by(location, Artifact_type) %>%
    mutate(percent = count/sum(count))
  return(df1)
}

ro_at = recycled_objects_artifact_type(artifacts)
# ggplot(ro_at) +
#   geom_col(aes(x = Artifact_type, y = percent, fill = recycled)) +
#   facet_wrap(~location) +
#   coord_flip()


##chi squared tests -- artifact type and recycling
test = table(artifacts %>% dplyr::select(location, Artifact_type, recycled))
dim(test)
dimnames(test)
test2 = apply(test, c(1,2,3), sum)

M = loglm( ~ Artifact_type * recycled * location, dat=test2, fit=TRUE)
M2 = update(M, ~ . - Artifact_type:location:recycled)
M3 = update(M, ~ . - Artifact_type:location:recycled - Artifact_type:location )
M4 = update(M, ~ . - Artifact_type:location:recycled - recycled:location)
M5 = update(M, ~ . - Artifact_type:location:recycled - recycled:Artifact_type)

anova(M, M2)
anova(M, M3)
anova(M, M4)
anova(M, M5)
#all the two way interactions seem to all contribute signficantly
#to the model, as does the three way interaction term -- interpretation?


#Raw materials by location
ggplot(artifacts) +
  geom_bar(aes(Raw_material_description)) +
  facet_wrap(~location) +
  coord_flip()

rml.table = table(artifacts %>% dplyr::select(location,Raw_material_description))
dim(rml.table)
dimnames(rml.table)
rml.table2 = apply(rml.table, c(1,2), sum)

M7 = loglm( ~ Raw_material_description * location, dat=rml.table2, fit=TRUE)
M8 = update(M7, ~ . - location:Raw_material_description)
anova(M7, M8)
chisq.test(rml.table2)
fisher.test(rml.table2, simulate.p.value = T)
#results indicate that there is a significant interaction between location and raw material type
#this means that there is a difference in the types of raw material amounts between locations
#this is consistent with visual findings, because P5 had many more different types of raw materials than P1 or P2

##chi squared tests -- Raw material description and recycling
test = table(artifacts %>% dplyr::select(location, Raw_material_description, recycled))
dim(test)
dimnames(test)
test2 = apply(test, c(1,2,3), sum)

M = loglm( ~ Raw_material_description * recycled * location, dat=test2, fit=TRUE)
M2 = update(M, ~ . - Raw_material_description:location:recycled)
M3 = update(M, ~ . - Raw_material_description:location:recycled - Raw_material_description:location )
M4 = update(M, ~ . - Raw_material_description:location:recycled - recycled:location)
M5 = update(M, ~ . - Raw_material_description:location:recycled - recycled:Raw_material_description)

anova(M, M2)
#three way interaction term does not contribute significantly to the model
anova(M, M3)
#raw material:location contributes significantly
anova(M, M4)
#recycled:location contributes signficantly
anova(M, M5)
#Raw_material:recycling interaction does not contribute significantly to the model

#within each site, raw material does not dictate whether or not something will be recycled


#Weathering class by location
ggplot(artifacts) +
  geom_bar(aes(Weathering_class)) +
  facet_wrap(~location) +
  coord_flip()

wcl.table = table(artifacts %>% dplyr::select(location,Weathering_class))
dim(wcl.table)
dimnames(wcl.table)
wcl.table2 = apply(wcl.table, c(1,2), sum)

M9 = loglm( ~ Weathering_class * location, dat=wcl.table2, fit=TRUE)
M10 = update(M9, ~ . - location:Weathering_class)
#anova(M9, M10)
chisq.test(wcl.table2)
#results indicate that there is a significant relationship between location and weathering class
#this means that there is a difference in the weathering class amounts between locations

##chi squared tests -- weathering class and recycling
test = table(artifacts %>% dplyr::select(location, Weathering_class, recycled))
dim(test)
dimnames(test)
test2 = apply(test, c(1,2,3), sum)

M = loglm( ~ Weathering_class * recycled * location, dat=test2, fit=TRUE)
M2 = update(M, ~ . - Weathering_class:location:recycled)
M3 = update(M, ~ . - Weathering_class:location )
M4 = update(M, ~ . - recycled:location)
M5 = update(M, ~ . - recycled:Weathering_class)

anova(M, M2)
#three way interaction contributes significantly to the model
anova(M, M3)
anova(M, M4)
anova(M, M5)
#so do all of the two way interactions
#this indicates that weathering class different affects which artifacts are recycled at different locations

#Size of recycled objects by location
plot(ks_test_recycled_vs_not(p1, positions = c(160, 175, 90, 2500)))
plot(ks_test_recycled_vs_not(p2, positions = c(150, 155, 80, 1050)))
plot(ks_test_recycled_vs_not(p5, positions = c(180, 160, 100, 1300)))

plot(ks_test_recycled_vs_not(artifacts))


dist_comp_recycled_vs_all(artifacts)


#Size of recycled objects compared between locations
dist_comp_recycled_vs_all(p1)
dist_comp_recycled_vs_all(p2)
dist_comp_recycled_vs_all(p5)

cowplot::plot_grid(dist_comp_recycled_vs_all(p1),
                   dist_comp_recycled_vs_all(p2),
                   dist_comp_recycled_vs_all(p5),
                   ncol = 1)

