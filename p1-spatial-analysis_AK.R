#Semizbugu P1 Analysis

library(tidyverse)
library(MASS)
library(rstatix)
library(ggpubr)
library(ggthemes)
library(gapminder)
library(ggridges)
library(cowplot)
library(lmtest)
library(lme4)
library(vcd)
library(ggstatsplot)
library(rcompanion)
library(DescTools)
library(writexl)

#point pattern data
library(sp)
library(rgdal)
library(spdep)
library(tidyverse)
library(sf)

#analysis of point pattern data
library(spatstat)
library(splancs)
library(maptools)


source("site-comparison-functions.R")
source("bordian_types_dictionary.R")

theme_set(theme_bw())

artifacts1 = read_csv("data/cleaned_june_artifacts.csv")
artifacts2 = read_csv("data/cleaned_july_artifacts.csv")
artifacts = rbind(artifacts1, artifacts2)
artifacts = artifacts %>% mutate(location = ifelse(Site_name %in% c("Square 1", "Square 2", "Square 3"), "Semizbugu P1", 
                                                   ifelse(Site_name %in% c("Square 4", "Square 5"), "Semizbugu P2", "Semizbugu P5")))
artifacts$recycled = !is.na(artifacts$Recycling_description)
artifacts$double_patina = str_detect(artifacts$Recycling_indications, "double_patina")
artifacts$Raw_material_description = tolower(artifacts$Raw_material_description)
artifacts = artifacts %>% filter(is.na(Problem_notes))

p1 = artifacts %>% filter(location == "Semizbugu P2")
p1$Weathering_class = factor(p1$Weathering_class, 
                                        levels = c("strongly_weathered", "mildly_weathered", "weakly_weathered", "not_weathered", "other"))
p1$Artifact_type = factor(p1$Artifact_type, 
                             levels = c("complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment", "shatter"))


sp.p1 = readOGR("data/artifact-shapefiles", "p2-artifacts")
sp.p1 = sp.p1[sp.p1$Id_nm %in% p1$Id_number,]


pal1 = c("#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
         "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

#summary stats
at.tab = as.data.frame(table(p1$Artifact_type)) %>%
  mutate(total = sum(Freq), 
         perc = Freq/total)


#size comparison
dist_comp_only_recycled(p1)
#ggsave(plot = dist_comp_only_recycled(p1), filename = "~/Desktop/NYU/Dissertation-Research/papers/Coco_AK/figures/p1-size-distributions-recycled.tiff")

p1.size = p1 %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 175 & Width <= 175 & Thickness <= 175) #size of calipers

ks.test((p1.size %>% filter(recycled == T))$Weight, p1.size$Weight, simulate.p.value = T)
ks.test((p1.size %>% filter(recycled == T))$Length, p1.size$Length, simulate.p.value = T)
ks.test((p1.size %>% filter(recycled == T))$Width, p1.size$Width, simulate.p.value = T)
ks.test((p1.size %>% filter(recycled == T))$Thickness, p1.size$Thickness, simulate.p.value = T)

#weathering comparison
ggplot(p1) +
  geom_bar(aes(x = Weathering_class, fill = recycled)) +
  scale_fill_colorblind()

p1.wc.r = as.data.frame(table(p1$recycled, p1$Weathering_class)) %>%
  group_by(Var2) %>%
  mutate(total = sum(Freq)) %>%
  mutate(perc = Freq/total) %>%
  mutate(perc = ifelse(perc == 0, NA, perc)) %>%
  mutate(Var2 = str_replace(Var2, "_", " ")) 
p1.wc.r$Var2 = factor(p1.wc.r$Var2, 
                             levels = c("strongly weathered", "mildly weathered", "weakly weathered", "not weathered", "other"))

wc.p1.plot = ggplot(p1.wc.r, aes(x=Var2, y=Freq, fill = Var1)) +
  geom_bar(stat = "identity", width=0.9) +
  geom_text(aes(label = paste0(round(perc*100),"%")), position = position_stack(vjust = 0.65), size = 2) +
  labs(x = "Weathering stage", y = "count", fill = "Recycled?") +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6")) 
plot(wc.p1.plot)
# ggsave(plot = wc.p1.plot, filename = "~/Desktop/NYU/Dissertation-Research/papers/Coco_AK/figures/p1-weathering-stage-recycled.tiff",
#        width = 7, height = 4.5)

#artifact type comparison
ggplot(p1) +
  geom_bar(aes(x = Artifact_type, fill = recycled)) +
  scale_fill_colorblind() +
  coord_flip()

p1.at.r = as.data.frame(table(p1$recycled, p1$Artifact_type)) %>%
  group_by(Var2) %>%
  mutate(total = sum(Freq)) %>%
  mutate(perc = Freq/total) %>%
  mutate(perc = ifelse(perc == 0, NA, perc)) %>%
  mutate(Var2 = str_replace(Var2, "_", " "))
p1.at.r$Var2 = factor(p1.at.r$Var2, 
                          levels = c("complete flake", "broken flake", "tool", "tool fragment", "core", "core fragment", "shatter"))


at.p1.plot = ggplot(p1.at.r %>% filter(!is.na(perc)), aes(x=Var2, y=Freq, fill = Var1)) +
  geom_bar(stat = "identity", width=0.9) +
  geom_text(aes(label = paste0(round(perc*100),"%")), position = position_stack(vjust = 0.70), size = 2) +
  coord_flip() +
  labs(x = "Artifact type", y = "count", fill = "Recycled?") +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"))
# ggsave(plot = at.p1.plot, filename = "~/Desktop/NYU/Dissertation-Research/papers/Coco_AK/figures/p1-artifact-type-recycled.tiff", 
#        width = 6, height = 3)


#tool type comparison
#need to do some editting of this -- make column for Bordian types
p1$Bordian_name = ""
for(i in 1:nrow(p1)) {
  print(p1$Bordian_type[i])
  if(!is.na(p1$Bordian_type[i])) {
    p1$Bordian_name[i] = bordian_types[p1$Bordian_type[i]]
  }
}

ggplot(p1 %>% filter(!is.na(Tool_type))) +
  geom_bar(aes(x = Bordian_name, fill = recycled)) +
  scale_fill_colorblind() +
  coord_flip()

p1.bt.r = as.data.frame(table(p1$recycled, p1$Bordian_name)) %>%
  group_by(Var2) %>%
  mutate(total = sum(Freq)) %>%
  mutate(perc = Freq/total) %>%
  mutate(perc = ifelse(perc == 0, NA, perc))

plot.bt = p1.bt.r %>% filter(!is.na(perc)) %>% filter(Var2 != "")
plot.bt$Var2 = factor(plot.bt$Var2, levels = rev(bordian_levels))

bt.p1.plot = ggplot(plot.bt, aes(x=Var2, y=Freq, fill = Var1)) +
  geom_bar(stat = "identity", width=0.9) +
  geom_text(aes(label = paste0(round(perc*100),"%")), position = position_stack(vjust = 0.65), size = 2) +
  coord_flip() +
  labs(x = "Bordian type", y = "count", fill = "Recycled?") +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6")) +
  theme(axis.text.y = element_text(size = 5))   
plot(bt.p1.plot)
# ggsave(plot = bt.p1.plot, filename = "~/Desktop/NYU/Dissertation-Research/papers/Coco_AK/figures/p1-bordian-type-recycled.tiff",
#        width = 7, height = 5)

#tool types
p1$Tool_type2 = ifelse(str_detect(p1$Tool_type, "notch denticulate"), "notch/denticulate", 
                       ifelse(str_detect(p1$Tool_type, " "), "multiple", 
                              p1$Tool_type))
unique(p1$Tool_type2)

p1.tt.r = as.data.frame(table(p1$recycled, p1$Tool_type2)) %>%
  group_by(Var2) %>%
  mutate(total = sum(Freq)) %>%
  mutate(perc = Freq/total) %>%
  mutate(perc = ifelse(perc == 0, NA, perc))

plot.tt = p1.tt.r %>% filter(!is.na(perc)) %>% filter(Var2 != "")
plot.tt$Var2 = factor(plot.tt$Var2, levels = rev(c("scraper", "notch", "denticulate", "notch/denticulate", 
                                               "biface", "point", "multiple", "other")))

tt.p1.plot = ggplot(plot.tt, aes(x=Var2, y=Freq, fill = Var1)) +
  geom_bar(stat = "identity", width=0.9) +
  geom_text(aes(label = paste0(round(perc*100),"%")), position = position_stack(vjust = 0.65), size = 2) +
  coord_flip() +
  labs(x = "Tool category", y = "count", fill = "Recycled?") +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"))
plot(tt.p1.plot)
# ggsave(plot = tt.p1.plot, filename = "~/Desktop/NYU/Dissertation-Research/papers/Coco_AK/figures/p1-tool-cat-recycled.tiff",
#        width = 6, height = 3.5)


# ggsave(plot = cowplot::plot_grid(at.p1.plot, tt.p1.plot, bt.p1.plot, ncol = 1, labels = "auto"), 
#        filename = "~/Desktop/NYU/Dissertation-Research/papers/Coco_AK/figures/p1-artifact-cat-recycled.tiff", 
#        width = 7.5, height = 10
# )


#recycling types
p1 = p1 %>%
  mutate(Recycling_type =
           ifelse(str_detect(Recycling_indications, "none"), "none", 
                  ifelse(str_detect(Recycling_indications, " "), "multiple", 
                                    Recycling_indications))) %>%
  mutate(Recycling_type = str_replace_all(Recycling_type, "_", " ")) %>%
  mutate(Recycling_type = ifelse(Recycling_type == "core on flake blade", "core on flake/blade", Recycling_type))
p1$Recycling_type = factor(p1$Recycling_type, 
                              levels = c("double patina", "core on flake/blade", "core on hammerstone", "multiple", "other", "none"))
unique(p1$Recycling_type)

ggplot(p1 %>% filter(!is.na(Recycling_type) & Recycling_type != "none")) +
  geom_bar(aes(x = Recycling_type)) +
  coord_flip()


#recycling types spatial
sp.p1@data = sp.p1@data %>%
  mutate(Recycling_type =
           ifelse(str_detect(Rcyclng_n, "none"), "none", 
                  ifelse(str_detect(Rcyclng_n, " "), "multiple", 
                         Rcyclng_n))) %>%
  mutate(Recycling_type = str_replace_all(Recycling_type, "_", " ")) %>%
  mutate(Recycling_type = ifelse(Recycling_type == "core on flake blade", "core on flake/blade", Recycling_type))
sp.p1@data$Recycling_type = factor(sp.p1@data$Recycling_type, 
                                   levels = c("double patina", "core on flake/blade", "core on hammerstone", "multiple", "other", "none"))

sf.p1 = st_as_sf(sp.p1)

pal1 = c("#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
         "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

plot.sf.p1 = sf.p1[!is.na(sf.p1$Recycling_type),]

rt.spat.plot = ggplot() +
  geom_sf(data = plot.sf.p1[plot.sf.p1$Recycling_type == "none",], size = 0.25, alpha = 0.25) +
  geom_sf(data = plot.sf.p1[plot.sf.p1$Recycling_type != "none",], aes(color = Recycling_type), size = 1) +
  #scale_color_colorblind() +
  scale_color_manual(values = pal1) +
  labs(color = "Recycling signature")
#ggsave(plot = rt.spat.plot, file = "~/Desktop/NYU/Dissertation-Research/papers/Coco_AK/figures/p1-recycling-signature-spatial.tiff")

plot(rt.spat.plot)

#regressions
#Artifact type
reg.data = p1 %>% select(recycled, Artifact_type)
logmodel = glm(recycled ~ Artifact_type, data = reg.data, family = binomial())
summary(logmodel)

reg.data = p1 %>% select(recycled, Bordian_name)  
logmodel = glm(recycled ~ Bordian_name, data = reg.data, family = binomial())
summary(logmodel)
unique(reg.data$Bordian_name)

reg.data = p1 %>% select(recycled, Tool_type2) %>%
  filter(!is.na(Tool_type2))
logmodel = glm(recycled ~ Tool_type2, data = reg.data, family = binomial())
summary(logmodel)  

reg.data = p1 %>% select(recycled, Weathering_class)
logmodel = glm(recycled ~ Weathering_class, data = reg.data, family = binomial())
summary(logmodel) 
exp(-0.5372)
exp(-1.4902)
exp(-1.6954)
