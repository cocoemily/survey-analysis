library(tidyverse)
library(ggpubr)
library(ggthemes)
library(MASS)

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


####assemblage history####

weibull.df = data.frame(
  location = character(), 
  w.shape = numeric(), 
  w.scale = numeric()
)
#####weibull distributions#####
p1.cf = all_artifacts %>% filter(location == "Semizbugu P1") %>%
  filter(Artifact_type == "complete_flake")
p1.pf = all_artifacts %>% filter(location == "Semizbugu P1") %>%
  filter(Artifact_type == "broken_flake" & Flake_fragment == "proximal")
p1.t = all_artifacts %>% filter(location == "Semizbugu P1") %>%
  filter(Artifact_type == "tool")
p1.tf = all_artifacts %>% filter(location == "Semizbugu P1") %>%
  filter(Artifact_type == "tool_fragment" & Flake_fragment == "proximal")
p1.w = rbind(p1.cf, p1.pf, p1.t, p1.tf)
p1.lengths = p1.w %>% filter(Length >= 20)
p1.lengths$new_Length = p1.lengths$Length - 19.99
hist(p1.lengths$new_Length)
w = fitdistr(p1.lengths$new_Length, densfun = "weibull")
weibull.df[nrow(weibull.df) + 1, ] = c("Semizbugu P1", w[["estimate"]]["shape"], w[["estimate"]]["scale"])

p2.cf = all_artifacts %>% filter(location == "Semizbugu P2") %>%
  filter(Artifact_type == "complete_flake")
p2.pf = all_artifacts %>% filter(location == "Semizbugu P2") %>%
  filter(Artifact_type == "broken_flake" & Flake_fragment == "proximal")
p2.t = all_artifacts %>% filter(location == "Semizbugu P2") %>%
  filter(Artifact_type == "tool")
p2.tf = all_artifacts %>% filter(location == "Semizbugu P2") %>%
  filter(Artifact_type == "tool_fragment" & Flake_fragment == "proximal")
p2.w = rbind(p2.cf, p2.pf, p2.t, p2.tf)
p2.lengths = p2.w %>% filter(Length >= 20)
p2.lengths$new_Length = p2.lengths$Length - 19.99
hist(p2.lengths$new_Length)
w = fitdistr(p2.lengths$new_Length, densfun = "weibull")
weibull.df[nrow(weibull.df) + 1, ] = c("Semizbugu P2", w[["estimate"]]["shape"], w[["estimate"]]["scale"])

p5.cf = all_artifacts %>% filter(location == "Semizbugu P5") %>%
  filter(Artifact_type == "complete_flake")
p5.pf = all_artifacts %>% filter(location == "Semizbugu P5") %>%
  filter(Artifact_type == "broken_flake" & Flake_fragment == "proximal")
p5.t = all_artifacts %>% filter(location == "Semizbugu P5") %>%
  filter(Artifact_type == "tool")
p5.tf = all_artifacts %>% filter(location == "Semizbugu P5") %>%
  filter(Artifact_type == "tool_fragment" & Flake_fragment == "proximal")
p5.w = rbind(p5.cf, p5.pf, p5.t, p5.tf)
p5.lengths = p5.w %>% filter(Length >= 20)
p5.lengths$new_Length = p5.lengths$Length - 19.99
hist(p5.lengths$new_Length)
w = fitdistr(p5.lengths$new_Length, densfun = "weibull")
weibull.df[nrow(weibull.df) + 1, ] = c("Semizbugu P5", w[["estimate"]]["shape"], w[["estimate"]]["scale"])

s4.cf = all_artifacts %>% filter(location == "Semizbugu 4") %>%
  filter(Artifact_type == "complete_flake")
s4.pf = all_artifacts %>% filter(location == "Semizbugu 4") %>%
  filter(Artifact_type == "broken_flake" & Flake_fragment == "proximal")
s4.t = all_artifacts %>% filter(location == "Semizbugu 4") %>%
  filter(Artifact_type == "tool")
s4.tf = all_artifacts %>% filter(location == "Semizbugu 4") %>%
  filter(Artifact_type == "tool_fragment" & Flake_fragment == "proximal")
s4.w = rbind(s4.cf, s4.pf, s4.t, s4.tf)
s4.lengths = s4.w %>% filter(Length >= 20)
s4.lengths$new_Length = s4.lengths$Length - 19.99
hist(s4.lengths$new_Length)
w = fitdistr(s4.lengths$new_Length, densfun = "weibull")
weibull.df[nrow(weibull.df) + 1, ] = c("Semizbugu 4", w[["estimate"]]["shape"], w[["estimate"]]["scale"])

s10a.cf = all_artifacts %>% filter(location == "Semizbugu 10A") %>%
  filter(Artifact_type == "complete_flake")
s10a.pf = all_artifacts %>% filter(location == "Semizbugu 10A") %>%
  filter(Artifact_type == "broken_flake" & Flake_fragment == "proximal")
s10a.t = all_artifacts %>% filter(location == "Semizbugu 10A") %>%
  filter(Artifact_type == "tool")
s10a.tf = all_artifacts %>% filter(location == "Semizbugu 10A") %>%
  filter(Artifact_type == "tool_fragment" & Flake_fragment == "proximal")
s10a.w = rbind(s10a.cf, s10a.pf, s10a.t, s10a.tf)
s10a.lengths = s10a.w %>% filter(Length >= 20)
s10a.lengths$new_Length = s10a.lengths$Length - 19.99
hist(s10a.lengths$new_Length)
w = fitdistr(s10a.lengths$new_Length, densfun = "weibull")
weibull.df[nrow(weibull.df) + 1, ] = c("Semizbugu 10A", w[["estimate"]]["shape"], w[["estimate"]]["scale"])

weibull.df$w.shape = as.numeric(weibull.df$w.shape)
weibull.df$w.scale = as.numeric(weibull.df$w.scale)

#### experimental weibull parameters ####
exp = read_csv("data/Experimental_chert_platform_flakes.csv")

ggplot(exp) +
  geom_histogram(aes(x = MAXLENGTH)) +
  facet_wrap(~Core.ID)

exp.weibull = exp %>% group_by(Core.ID) %>%
  filter(MAXLENGTH >= 20) %>%
  mutate(newlength = MAXLENGTH - 19.99) %>%
  summarize(
    w.shape = fitdistr(newlength, densfun = "weibull")[["estimate"]]["shape"], 
    w.scale = fitdistr(newlength, densfun = "weibull")[["estimate"]]["scale"]
  )

ggplot(exp.weibull) +
  geom_point(aes(x = w.scale, y = w.shape))
  
#### archaeological weibull parameters ####
pech = data.frame(
  location = c("Layer 3A", "Layer 3B", "Layer 4", "Layer 5A", "Layer 5B", "Layer 6A", "Layer 6B", "Layer 7", "Layer 8"),
  w.shape = c(1.16, 1.14, 1.14, 1.20, 1.19, 1.21, 1.21, 1.16, 1.22), 
  w.scale = c(10.61, 10.59, 15.10, 14.18, 14.63, 16.07, 10.72, 10.86, 13.07)
)


wplot = ggplot(mapping = aes(x = w.scale, y = w.shape)) +
  geom_point(data = exp.weibull, color = "grey", mapping = aes(shape = "experimental")) +
  geom_point(data = pech, mapping = aes(shape = "Pech IV"), size = 3) +
  #geom_point(data = aotearoa, mapping = aes(shape = "Aotearoa"), size = 3) +
  geom_point(data = weibull.df, mapping = aes(color = location), size = 3) +
  scale_color_tableau() +
  labs(x = "Weibull scale (distribution spread)", 
       y = "Weibull shape (distribution slope)") +
  scale_shape_manual(name = "",
                     breaks = c("experimental", "Pech IV", "Semizbugu"), 
                     values = c("experimental" = 15, "Pech IV" = 3, "Semizbugu" = 16))

ggsave(filename = "figures/weibull-comparison.tiff", wplot, 
       dpi = 300, width = 7, height = 5)
  