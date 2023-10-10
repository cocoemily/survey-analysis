##flake size analysis
library(tidyverse)
library(ggpubr)
library(ggthemes)

set.seed(11122)

theme_set(theme_bw())

#### DATA LOADING ####
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
  filter(is.na(Length) | Length <= 200) %>%
  filter(is.na(Width) | Width <= 200) %>% 
  filter(is.na(Thickness) |Thickness <= 200) #size of calipers

all_artifacts = all_artifacts %>%
  filter(Weight <= 1000)

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


#### platform area calculations ####
complete_flakes = all_artifacts %>% filter(Artifact_type == "complete_flake") %>% 
  filter(Retouch == FALSE) %>% filter(recycled == FALSE)
complete_flakes$platform_area = complete_flakes$Platform_thickness * complete_flakes$Platform_width
complete_flakes$surface_area = complete_flakes$Length * complete_flakes$Width

ggplot(complete_flakes) +
  geom_point(aes(x = log(platform_area), y = log(Weight))) +
  geom_smooth(aes(x = log(platform_area), y = log(Weight)), method = "lm") +
  scale_color_colorblind()

fit1 = lm(log(Weight) ~ log(platform_area), data = complete_flakes)
summary(fit1)

fit2 = lm(log(surface_area) ~ log(platform_area), data = complete_flakes)
summary(fit2)


#### predict flake size ####
all_cflakes = all_artifacts %>%
  filter(!is.na(Platform_thickness) & !is.na(Platform_width))
all_cflakes$platform_area = all_cflakes$Platform_thickness * all_cflakes$Platform_width

all_cflakes$predicted_weight = exp(predict(fit1, newdata = all_cflakes))
all_cflakes$predicted_sa = exp(predict(fit2, newdata = all_cflakes))

predicted = all_cflakes %>% pivot_longer(
  cols = c(predicted_weight, predicted_sa), 
  names_to = "attribute", values_to = "predicted_val"
)
predicted$attribute = factor(predicted$attribute, 
                             levels = c("predicted_weight", "predicted_sa"))

p.labs = c("weight", "surface area")
names(p.labs) = c("predicted_weight", "predicted_sa")

weights.r = (predicted %>% filter(recycled == T) %>% filter(attribute == "predicted_weight"))$predicted_val
weights.u = (predicted %>% filter(recycled == F) %>% filter(attribute == "predicted_weight"))$predicted_val

hist(weights.r)
hist(weights.u)

pred_w = predicted %>% filter(attribute == "predicted_weight")
wilcox.test(predicted_val ~ recycled, data = pred_w, alternative = "less")

pred_sa = predicted %>% filter(attribute == "predicted_sa")
wilcox.test(predicted_val ~ recycled, data = pred_sa, alternative = "less")


p1 = ggplot(predicted, aes(y = predicted_val, x = as.factor(recycled), 
                        color = as.factor(recycled))) +
  geom_jitter(alpha = 0.2, size = 0.5) +
  geom_boxplot(alpha = 0.4) +
  facet_wrap(~attribute, scales = "free_y", labeller = labeller(attribute=p.labs)) +
  stat_compare_means(label = "p.signif", label.x.npc = "middle", vjust = 2) +
  scale_x_discrete(labels = c("unrecycled", "recycled")) +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  labs(x = "", y = "predicted value", color = "") +
  guides(color = "none") +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7))
plot(p1)
ggsave(filename = "figures/predicted_flake_sizes.tiff", p1, 
       dpi = 300, width = 6, height = 4)

p2 = ggplot(predicted %>% filter(location %in% c("Semizbugu P1", "Semizbugu P2", "Semizbugu P5")), 
            aes(y = predicted_val, x = as.factor(recycled), 
                           color = as.factor(recycled))) +
  geom_jitter(alpha = 0.2, size = 0.5) +
  geom_boxplot(alpha = 0.4) +
  facet_wrap(~attribute, scales = "free_y", labeller = labeller(attribute=p.labs)) +
  stat_compare_means(label = "p.signif", label.x.npc = "middle", vjust = 2) +
  scale_x_discrete(labels = c("unrecycled", "recycled")) +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  labs(x = "", y = "predicted value", color = "") +
  guides(color = "none") +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7))
plot(p2)

p.supp = ggplot(predicted, aes(y = predicted_val, x = as.factor(recycled), 
                      color = as.factor(recycled))) +
  geom_jitter(alpha = 0.2, size = 0.5) +
  geom_boxplot(alpha = 0.4) +
  facet_grid(attribute~location, scales = "free", labeller = labeller(attribute=p.labs)) +
  stat_compare_means(label = "p.signif", label.x.npc = "middle", vjust = 2) +
  scale_x_discrete(labels = c("unrecycled", "recycled")) +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  labs(x = "", y = "predicted value", color = "") +
  guides(color = "none") +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7))
plot(p.supp)
ggsave(filename = "figures/predicted_flake_sizes_by-location.tiff", p.supp, 
       dpi = 300, width = 9, height = 4)

p.recycled = all_cflakes %>% filter(recycled == TRUE)
p.nrecycled = all_cflakes %>% filter(recycled == FALSE)

wilcox.test(p.recycled$predicted_weight, p.nrecycled$predicted_weight, alternative = "greater")
wilcox.test(p.recycled$predicted_sa, p.nrecycled$predicted_sa, alternative = "greater")


rcompanion::plotNormalHistogram(sqrt(p.recycled$predicted_weight))
rcompanion::plotNormalHistogram(sqrt(p.recycled$predicted_sa))
