#comparison of different sites
library(tidyverse)
library(MASS)
library(rstatix)
library(ggpubr)
library(ggthemes)
library(gapminder)
library(ggridges)
library(lmtest)
library(vcd)
library(ggstatsplot)
library(rcompanion)
library(lme4)
library(QuantPsyc)
library(broom.mixed)

set.seed(11122)

source("site-comparison-functions.R")

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


#### TO DO: raw material cleaning ####

#possible artifact rolling
all_artifacts$poss_roll = ifelse(all_artifacts$Bordian_name == "alternate scraper", TRUE, FALSE)
all_artifacts$poss_roll = ifelse(str_detect(all_artifacts$Retouch_side, pattern = " "), TRUE, all_artifacts$poss_roll)

rolled.rcycl = all_artifacts %>% filter(poss_roll == T & recycled == T)
all_artifacts = subset(all_artifacts, !(Id_number %in% rolled.rcycl$Id_number))


#### cleaning up environment ####
rm(list = c("artifacts", "artifacts1", "artifacts2", "collections", "paleocore_sss_artifact_form_all_versions_False_2022_08_01_05_36_18", "s10", "rolled.rcycl"))


####recycling intensity values####
all.p1 = nrow(all_artifacts %>% filter(location == "Semizbugu P1"))
recycl.p1 = nrow(all_artifacts %>% filter(location == "Semizbugu P1") %>% filter(recycled == T))
ri.p1 = recycl.p1/all.p1

all.p2 = nrow(all_artifacts %>% filter(location == "Semizbugu P2"))
recycl.p2 = nrow(all_artifacts %>% filter(location == "Semizbugu P2") %>% filter(recycled == T))
ri.p2 = recycl.p2/all.p2

all.p5 = nrow(all_artifacts %>% filter(location == "Semizbugu P5"))
recycl.p5 = nrow(all_artifacts %>% filter(location == "Semizbugu P5") %>% filter(recycled == T))
ri.p5 = recycl.p5/all.p5

all.s4 = nrow(all_artifacts %>% filter(location == "Semizbugu 4"))
recycl.s4 = nrow(all_artifacts %>% filter(location == "Semizbugu 4") %>% filter(recycled == T))
ri.s4 = recycl.s4/all.s4

all.s10a = nrow(all_artifacts %>% filter(location == "Semizbugu 10A"))
recycl.s10a = nrow(all_artifacts %>% filter(location == "Semizbugu 10A") %>% filter(recycled == T))
ri.s10a = recycl.s10a/all.s10a

####recycling intensity correlations####
ri.cor = data.frame(
  location = unique(all_artifacts$location),
  ri = c(ri.p1, ri.p2, ri.p5, ri.s10a, ri.s4),
  r.count = c(recycl.p1, recycl.p2, recycl.p5, recycl.s10a, recycl.s4),
  count = c(all.p1, all.p2, all.p5, all.s10a, all.s4), 
  cr = c(0.54, 0.65, 0.45, 0.46, 0.50)
)

cplot1 = ggplot(ri.cor) +
  #geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = log(count), y = ri, color = location)) +
  geom_smooth(aes(x = log(count), y = ri), method = "lm", se = T) +
  scale_color_tableau() +
  labs(x = "log(artifact count)", y = "recycling intensity") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75))
plot(cplot1)

cplot2 = ggplot(ri.cor) +
  #geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = cr, y = ri, color = location)) +
  geom_smooth(aes(x = cr, y = ri), method = "lm", se = T) +
  scale_color_tableau() +
  labs(x = "cortex ratio", y = "recycling intensity") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75))
plot(cplot2)

# ggsave(
#   filename = "figures/recycling-intensity-correlations.tiff",
#   ggarrange(cplot1, cplot2, ncol = 2, nrow = 1, common.legend = T, labels = "AUTO", legend = "bottom"), 
#   dpi = 300, width = 8, height = 4
# )

#### Differences by location ####
ggplot(all_artifacts) +
  geom_bar(aes(recycled, fill = recycled)) +
  facet_wrap(~location, ncol = 1) +
  coord_flip() +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6")) +
  theme(legend.position = "bottom")


rl.table = table(all_artifacts %>% dplyr::select(location,recycled))

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
at.r.df = as.data.frame(Table)
at.r.cmh = groupwiseCMH(Table,
                        group   = 1,
                        fisher  = TRUE,
                        gtest   = FALSE,
                        chisq   = FALSE,
                        method  = "fdr",
                        correct = "none",
                        digits  = 3, 
                        simulate.p.value = T)
at.r.df = at.r.df %>% left_join(at.r.cmh, by = c("location" = "Group"))
# Hypotheses
# •  Null hypothesis:  There is no association between the two inner variables.
# •  Alternative hypothesis (two-sided): There is an association between the two inner variables.
##suggests that there an association between recycling and artifact type at all sites, except for Semizbugu 10A

at.r.df = at.r.df %>% group_by(location, Artifact_type) %>%
  mutate(loc.total = sum(Freq), 
         percentage = Freq/loc.total, 
         Artifact_type = str_replace(Artifact_type, "_", " "))
at.r.df$Artifact_type = factor(at.r.df$Artifact_type , 
                               levels = c("complete flake", "broken flake", "tool", "tool fragment", "core", "core fragment", "shatter"))
at.r.df$location = factor(at.r.df$location, levels = c("Semizbugu P1", "Semizbugu P2", "Semizbugu P5", "Semizbugu 10A", "Semizbugu 4"))
#####recycling by artifact type ######
at.r.plot = ggplot(at.r.df, aes(x = Artifact_type, y = Freq, fill = recycled)) +
  geom_bar(stat = "identity", width=0.9, position = "dodge") +
  geom_text(data = at.r.df %>% filter(percentage != 0), aes(label = paste0(round(percentage*100, digits = 2),"%")), position = position_dodge2(width = 1), size = 2, hjust = -0.01) +
  geom_text(aes(label = paste("P = ", adj.p), x = "shatter", y = 450), size = 3, fontface = "italic", hjust = "right", vjust = -0.1) +
  facet_wrap(~location, strip.position = "top") +
  coord_flip() +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7), 
        legend.position = "bottom") +
  labs(x = "", y = "count", fill = "")
plot(at.r.plot)

ggsave(filename = "figures/artifact-type-recycling.tiff", at.r.plot, 
       dpi = 300, width = 10, height = 6)

all_artifacts$Artifact_type = factor(all_artifacts$Artifact_type, 
                                     levels = c(
                                       "shatter", "complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment"
                                     ))
rafit = glmer(recycled ~ Artifact_type + (1 | location), family = binomial(), data = all_artifacts)
summary(rafit)

wc.table = table(all_artifacts %>% dplyr::select(location, Weathering_class) %>%
                   filter(!is.na(Weathering_class)))
pairwiseNominalIndependence(wc.table, simulate.p.value = T, 
                            fisher = T, chisq = F, gtest = F)
#no difference between P2 and P5 in distribution of artifacts among weathering classes

wc.r.table = table(all_artifacts %>% dplyr::select(location, recycled, Weathering_class) %>%
                     filter(!is.na(Weathering_class)))
wc.r.table2 = apply(wc.r.table, c(1,2,3), sum)
wc.r = as.data.frame(ftable(wc.r.table))
Table = xtabs(Freq ~ location + recycled + Weathering_class, data=wc.r)
wc.r.df = as.data.frame(Table)
wc.r.cmh = groupwiseCMH(Table,
                        group   = 1,
                        fisher  = TRUE,
                        gtest   = FALSE,
                        chisq   = FALSE,
                        method  = "fdr",
                        correct = "none",
                        digits  = 3, 
                        simulate.p.value = T)
wc.r.df = wc.r.df %>% left_join(wc.r.cmh, by = c("location" = "Group"))
# Hypotheses
# •  Null hypothesis:  There is no association between the two inner variables.
# •  Alternative hypothesis (two-sided): There is an association between the two inner variables.
##suggests that there an association between recycling and weathering class at all sites, except for Semizbugu P2

wc.r.df = wc.r.df %>% group_by(location, Weathering_class) %>%
  mutate(loc.total = sum(Freq), 
         percentage = Freq/loc.total, 
         Weathering_class = str_replace(Weathering_class, "_", " "))
wc.r.df$Weathering_class = factor(wc.r.df$Weathering_class , 
                                  levels = c("strongly weathered", "mildly weathered", "weakly weathered", "not weathered", "other"))
wc.r.df$location = factor(wc.r.df$location, levels = c("Semizbugu P1", "Semizbugu P2", "Semizbugu P5", "Semizbugu 10A", "Semizbugu 4"))
#####recycling by weathering class ######
wc.r.plot = ggplot(wc.r.df, aes(x = Weathering_class, y = Freq, fill = recycled)) +
  geom_bar(stat = "identity", width=0.9, position = "dodge") +
  geom_text(data = wc.r.df %>% filter(percentage != 0), aes(label = paste0(round(percentage*100, digits = 2),"%")), position = position_dodge2(width = 1), size = 2, hjust = -0.01) +
  geom_text(aes(label = paste("P = ", adj.p), x = "other", y = 500), size = 3, fontface = "italic", hjust = "right", vjust = -1.5) +
  facet_wrap(~location, strip.position = "top") +
  coord_flip() +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7), 
        legend.position = "bottom") +
  labs(x = "", y = "count", fill = "")
plot(wc.r.plot)

ggsave(filename = "figures/weathering-class-recycling.tiff", wc.r.plot, 
       dpi = 300, width = 10, height = 6)

wc.r.plot2 = ggplot(wc.r.df, aes(x = Weathering_class, y = Freq, fill = recycled)) +
  geom_bar(stat = "identity", width=0.9, position = "dodge") +
  geom_text(data = wc.r.df %>% filter(percentage != 0), aes(label = paste0(round(percentage*100, digits = 2),"%")), position = position_dodge2(width = 1), size = 3, hjust = -0.1) +
  geom_text(aes(label = paste("P = ", adj.p), x = "other", y = 500), size = 4.5, fontface = "italic", hjust = "right", vjust = -1.5) +
  facet_wrap(~location, nrow = 1, strip.position = "top") +
  coord_flip() +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6")) +
  theme(strip.text = element_text(size = 9, face = "bold"), axis.text = element_text(size = 8)) +
  labs(x = "weathering class", y = "count", fill = "recycled?")
plot(wc.r.plot2)
ggsave(filename = "figures/SAA_weathering-class-recycling.tiff", wc.r.plot2, 
       dpi = 300, width = 13, height = 6)

wreg = all_artifacts
wreg$Weathering_class = factor(wreg$Weathering_class, levels = c("not_weathered", "strongly_weathered", "mildly_weathered", "weakly_weathered", "other"))
rwfit = glmer(recycled ~ Weathering_class + (1 | location), family = binomial(), data = wreg)
rw.df = tidy(rwfit)

rw.df$term_clean = c("(intercept)", "strongly weathered", "mildly weathered", "weakly weathered", "other", "")
rw.df$term_clean = factor(rw.df$term_clean, levels = c("(intercept)", "strongly weathered", "mildly weathered", "weakly weathered", "other", ""))
rw.df$lower = rw.df$estimate - rw.df$std.error
rw.df$upper = rw.df$estimate + rw.df$std.error

lmerrw = ggplot(rw.df %>% filter(effect == "fixed")) +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  geom_point(aes(x = term_clean,  y = estimate), size = 6, color = "blue") +
  geom_linerange(aes(x = term_clean,  y = estimate, ymin = lower, ymax = upper), color = "blue") +
  coord_flip() +
  labs(x = "term")
ggsave(filename = "figures/SAA_weathering-class-recycling_lmer.tiff", lmerrw , 
       dpi = 300, width = 6, height = 5)


# pairwiseNominalIndependence(
#   table(all_artifacts %>% dplyr::select(location, Retouch)),
#   simulate.p.value = T, 
#   fisher = T, chisq = F, gtest = F
# )
#Semizbugu 4 has significant differences in distribution of retouch (presence/absence) compared to all other sites

# pairwiseNominalIndependence(
#   table(all_artifacts %>% dplyr::select(location, Edge_damage)),
#   simulate.p.value = T, 
#   fisher = T, chisq = F, gtest = F
# )

# pairwiseNominalIndependence(
#   table(all_artifacts %>% dplyr::select(location, retouch.side)),
#   simulate.p.value = T, 
#   fisher = T, chisq = F, gtest = F
# )

# rs.r.table = table(all_artifacts %>% dplyr::select(location, recycled, retouch.side))
# rs.r.table2 = apply(rs.r.table, c(1,2,3), sum)
# rs.r = as.data.frame(ftable(rs.r.table))
# Table = xtabs(Freq ~ location + recycled + retouch.side,
#               data=rs.r)
# groupwiseCMH(Table,
#              group   = 1,
#              fisher  = TRUE,
#              gtest   = FALSE,
#              chisq   = FALSE,
#              method  = "fdr",
#              correct = "none",
#              digits  = 3, 
#              simulate.p.value = T)
# # Hypotheses
# # •  Null hypothesis:  There is no association between the two inner variables.
# # •  Alternative hypothesis (two-sided): There is an association between the two inner variables.
# ##no association between recycling and retouch side at any sites

# pairwiseNominalIndependence(
#   table(all_artifacts %>% dplyr::select(location, Bordian_name)),
#   simulate.p.value = T, 
#   fisher = T, chisq = F, gtest = F
# )

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
tt.r.df = as.data.frame(Table)
tt.r.cmh = groupwiseCMH(Table,
                        group   = 1,
                        fisher  = TRUE,
                        gtest   = FALSE,
                        chisq   = FALSE,
                        method  = "fdr",
                        correct = "none",
                        digits  = 3, 
                        simulate.p.value = T)
tt.r.df = tt.r.df %>% left_join(tt.r.cmh, by = c("location" = "Group"))
# Hypotheses
# •  Null hypothesis:  There is no association between the two inner variables.
# •  Alternative hypothesis (two-sided): There is an association between the two inner variables.
##suggests that there an association between recycling and tool type at P5 and S4, but not at the other sites
tt.r.df = tt.r.df %>% group_by(location, tool.type) %>%
  mutate(loc.total = sum(Freq), 
         percentage = Freq/loc.total)
tt.r.df$tool.type = factor(tt.r.df$tool.type , 
                           levels = rev(c("scraper", "notch", "denticulate", "notch/denticulate", 
                                          "biface", "point", "multiple", "other")))
tt.r.df$location = factor(tt.r.df$location, levels = c("Semizbugu P1", "Semizbugu P2", "Semizbugu P5", "Semizbugu 10A", "Semizbugu 4"))
#####recycling by tool type ######
tt.r.plot = ggplot(tt.r.df, aes(x = tool.type, y = Freq, fill = recycled)) +
  geom_bar(stat = "identity", width=0.9, position = "dodge") +
  geom_text(data = tt.r.df %>% filter(percentage != 0), aes(label = paste0(round(percentage*100, digits = 2),"%")), position = position_dodge2(width = 1), size = 2, hjust = -0.01) +
  geom_text(aes(label = paste("P = ", adj.p), y = 175, x = "notch"), size = 3, fontface = "italic", hjust = "right") +
  facet_wrap(~location, strip.position = "top") +
  coord_flip() +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7), 
        legend.position = "bottom") +
  labs(x = "", y = "count", fill = "")
ggsave(filename = "figures/tool-type-recycling.tiff", tt.r.plot, 
       dpi = 300, width = 10, height = 6)

rtfit = glmer(recycled ~ tool.type + (1 | location), family = binomial(), data = all_artifacts %>% filter(!is.na(tool.type)))
summary(rtfit)
levels(all_artifacts$tool.type)

p5fit = glm(recycled ~ tool.type, family = binomial(), data = all_artifacts %>% filter(!is.na(tool.type)) %>% filter(location == "Semizbugu P5"))
summary(p5fit)
s4fit = glm(recycled ~ tool.type, family = binomial(), data = all_artifacts %>% filter(!is.na(tool.type)) %>% filter(location == "Semizbugu 4"))
summary(s4fit)

######p5 scrapers#####
p5.scrapers = all_artifacts %>% filter(location == "Semizbugu P5") %>%
  filter(tool.type == "scraper")
table(p5.scrapers$Weathering_class)


# pairwiseNominalIndependence(
#   table(all_artifacts %>% dplyr::select(location, flake.type)),
#   simulate.p.value = T, 
#   fisher = T, chisq = F, gtest = F
# )
# 
# pairwiseNominalIndependence(
#   table(all_artifacts %>% dplyr::select(location, Flake_fragment)),
#   simulate.p.value = T, 
#   fisher = T, chisq = F, gtest = F
# )
# 
# pairwiseNominalIndependence(
#   table(all_artifacts %>% dplyr::select(location, Flake_termination)),
#   simulate.p.value = T, 
#   fisher = T, chisq = F, gtest = F
# )

mfig1 = ggarrange(at.r.plot, tt.r.plot, wc.r.plot, nrow = 1, common.legend = T, legend = "bottom", labels = "AUTO")
ggsave(filename = "figures/cat-data-comparison-by-site.tiff", mfig1, 
       dpi = 300, width = 12, height = 8)


#####recycling by weight #####
pairwiseKS(all_artifacts %>% dplyr::select(location, Weight))
#statistically significant interaction between location and recycling on weight
pwt = all_artifacts %>% filter(Weight <= 1000) %>% group_by(location) %>%
  pairwise_wilcox_test(Weight ~ recycled, p.adjust.method = "bonferroni")
we.r.plot = ggplot(all_artifacts %>% filter(Weight <= 1000)) +
  geom_density(aes(Weight, group = recycled, fill = recycled, color = recycled), alpha = 0.5) +
  facet_wrap(~location, ncol = 1, strip.position = "right") +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  geom_text(data = pwt, aes(label = paste("P =",p.adj), x = 900, y = 0.025), fontface = "italic", hjust = "right") +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7)) +
  labs(x = "weight", y = "density", fill = "", color = "")
plot(we.r.plot)
ggsave(filename = "figures/weight-recycling.tiff", we.r.plot, 
       dpi = 300, width = 7, height = 8)

#####recycling by length #####
pairwiseKS(all_artifacts %>% dplyr::select(location, Length))
plt = all_artifacts  %>% group_by(location) %>%
  pairwise_wilcox_test(Length ~ recycled, p.adjust.method = "bonferroni")
l.r.plot = ggplot(all_artifacts) +
  geom_density(aes(Length, group = recycled, fill = recycled, color = recycled), alpha = 0.5) +
  facet_wrap(~location, ncol = 1, strip.position = "right") +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  geom_text(data = plt, aes(label = paste("P =",p.adj), x = 175, y = 0.020), fontface = "italic", hjust = "right") +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7)) +
  labs(x = "length", y = "density", fill = "", color = "")
ggsave(filename = "figures/length-recycling.tiff", l.r.plot, 
       dpi = 300, width = 7, height = 8)

#####recycling by width #####
pairwiseKS(all_artifacts %>% dplyr::select(location, Width))
pwdt = all_artifacts  %>% group_by(location) %>%
  pairwise_wilcox_test(Width ~ recycled, p.adjust.method = "bonferroni")
wd.r.plot = ggplot(all_artifacts) +
  geom_density(aes(Width, group = recycled, fill = recycled, color = recycled), alpha = 0.5) +
  facet_wrap(~location, ncol = 1, strip.position = "right") +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  geom_text(data = pwdt, aes(label = paste("P =",p.adj), x = 180, y = 0.025), fontface = "italic", hjust = "right") +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7)) +
  labs(x = "width", y = "density", fill = "", color = "")
ggsave(filename = "figures/width-recycling.tiff", wd.r.plot, 
       dpi = 300, width = 7, height = 8)


#####recycling by thickness #####
pairwiseKS(all_artifacts %>% dplyr::select(location, Thickness))
ptt = all_artifacts  %>% group_by(location) %>%
  pairwise_wilcox_test(Thickness ~ recycled, p.adjust.method = "bonferroni")
t.r.plot = ggplot(all_artifacts) +
  geom_density(aes(Thickness, group = recycled, fill = recycled, color = recycled), alpha = 0.5) +
  facet_wrap(~location, ncol = 1, strip.position = "right") +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  geom_text(data = ptt, aes(label = paste("P =",p.adj), x = 110, y = 0.1), fontface = "italic", hjust = "right") +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7)) +
  labs(x = "thickness", y = "density", fill = "", color = "")
ggsave(filename = "figures/thickness-recycling.tiff", t.r.plot, 
       dpi = 300, width = 7, height = 8)

mfig2 = ggarrange(l.r.plot, wd.r.plot, t.r.plot, we.r.plot, 
                  nrow = 1,
                  common.legend = T, legend = "bottom", labels = "AUTO")
ggsave(filename = "figures/cont-data-comparison-by-site1.tiff", mfig2, 
       dpi = 300, width = 13, height = 7)


#####recycling by cortex percentages#####
pairwiseKS(all_artifacts %>% dplyr::select(location, Cortex_percentage))
pcrt = all_artifacts  %>% group_by(location) %>%
  pairwise_wilcox_test(Cortex_percentage ~ recycled, p.adjust.method = "bonferroni")
cr.r.plot = ggplot(all_artifacts) +
  geom_density(aes(Cortex_percentage, group = recycled, fill = recycled, color = recycled), alpha = 0.5) +
  facet_wrap(~location, ncol = 1, strip.position = "right") +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  geom_text(data = pcrt, aes(label = paste("P =",p.adj), x = 95, y = 0.065), fontface = "italic", hjust = "right") +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7)) +
  labs(x = "cortex percentage", y = "density", fill = "", color = "")
ggsave(filename = "figures/cortex-ratio-recycling.tiff", cr.r.plot, 
       dpi = 300, width = 7, height = 8)

# pairwiseKS(all_artifacts %>% dplyr::select(location, Platform_thickness))
# pairwiseKS(all_artifacts %>% dplyr::select(location, Platform_width))

#####recycling by flake scar counts #####
pairwiseKS(all_artifacts %>% dplyr::select(location, Dorsal_flake_scar_count))
pdft = all_artifacts  %>% group_by(location) %>%
  pairwise_wilcox_test(Dorsal_flake_scar_count ~ recycled, p.adjust.method = "bonferroni")
df.r.plot = ggplot(all_artifacts) +
  geom_density(aes(Dorsal_flake_scar_count, group = recycled, fill = recycled, color = recycled), alpha = 0.5) +
  facet_wrap(~location, ncol = 1, strip.position = "right") +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("unrecycled", "recycled")) +
  geom_text(data = pdft, aes(label = paste("P =",p.adj), x = 55, y = 0.25), fontface = "italic", hjust = "right") +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 7)) +
  labs(x = "dorsal flake scar count", y = "density", fill = "", color = "")
ggsave(filename = "figures/flake-scar-count-recycling.tiff", df.r.plot, 
       dpi = 300, width = 7, height = 8)

mfig3 = ggarrange(cr.r.plot, df.r.plot, 
                  nrow = 1,
                  common.legend = T, legend = "bottom", labels = "AUTO")
ggsave(filename = "figures/cont-data-comparison-by-site2.tiff", mfig3, 
       dpi = 300, width = 7, height = 7)


#### Regressions between locations ####
##does location explain more variation than other variables?
reg.data = all_artifacts %>% dplyr::select("recycled","Weathering_class", 
                                           "Dorsal_flake_scar_count", "Cortex_percentage",
                                           "Thickness", "Length", "Width", 
                                           "Weight")
reg.data$recycled = as.factor(reg.data$recycled)
reg.data$Weathering_class = factor(reg.data$Weathering_class, levels = c("not_weathered", "strongly_weathered", "mildly_weathered", "weakly_weathered", "other"))

fit1 = glm(recycled ~ ., family = binomial(), data = reg.data)
summary(fit1)

reg.data2 = all_artifacts %>% dplyr::select("location", "recycled", "Weathering_class", 
                                            "Dorsal_flake_scar_count", "Cortex_percentage",
                                            "Thickness", "Length", "Width", 
                                            "Weight")
fit2 = glm(recycled ~ ., family = binomial(), data = reg.data2)
summary(fit2)

anova(fit1, fit2, test = "LR")
##need to include location, but weakly and not weathered artifacts are still less recycled, 
##and overall long and wide artifacts are more likely to be recycled, but also thinner artifacts

###### standardized coefficients#####
coef = as.data.frame(summary(fit2)$coefficients)
coef$var = rownames(coef)
scoef = as.data.frame(lm.beta::lm.beta(fit2)[["standardized.coefficients"]])
scoef$var = rownames(scoef)
colnames(scoef) = c("coef", "var")
rownames(scoef) = NULL
scoef$abs_coef = abs(scoef$coef)
scoef = scoef %>% left_join(coef)
scoef$signf = ifelse(scoef$`Pr(>|z|)` < 0.05, "TRUE", "FALSE")

scoef$var = factor(scoef$var, 
                   levels = c(
                     "locationSemizbugu 4", "locationSemizbugu P1", "locationSemizbugu P2", "locationSemizbugu P5", 
                     "Weathering_classmildly_weathered", "Weathering_classweakly_weathered", "Weathering_classnot_weathered", "Weathering_classother",
                     "Dorsal_flake_scar_count", "Cortex_percentage", 
                     "Length", "Width", "Thickness", "Weight"
                   ))
scoef$var = fct_rev(scoef$var)
sc.labs = rev(c("location: Semizbugu 4", "location: Semizbugu P1", "location: Semizbugu P2", "location: Semizbugu P5", 
                "weathering class: mild", "weathering class: weak", "weathering class: not", "weathering class: other",
                "dorsal flake scar count", "cortex percentage", 
                "length", "width", "thickness", "weight"))

sc.plot = ggplot(scoef %>% filter(!is.na(var))) +
  geom_hline(yintercept = 0, color = I("red")) +
  geom_hline(yintercept = min((scoef %>% filter(str_detect(var, "location")))$abs_coef), color = I("grey"), linetype = "dashed") +
  geom_point(aes(x = var, y = abs_coef, color = var, shape = signf), size = 5)+
  coord_flip() +
  guides(color = "none") +
  labs(x = "variable", y = "absolute value of standardized coefficient", 
       shape = "significant?") +
  scale_x_discrete(labels = sc.labs) +
  scale_shape_manual(values = c(4, 16))
plot(sc.plot)
ggsave(filename = "figures/stand-coefs-continuous-vars.tiff", sc.plot,
       dpi = 300, width = 8, height = 5)


#####multiple regressions for each of the locations #####
reg.data.p1 = all_artifacts %>% filter(location == "Semizbugu P1") %>%
  dplyr::select("recycled", "Weathering_class", "Artifact_type", "Tool_type", "Cortex_percentage",
                "Length", "Width", "Thickness", 
                "Weight")
reg.data.p1$Tool_type = as.character(reg.data.p1$Tool_type)
reg.data.p1$Tool_type = ifelse(is.na(reg.data.p1$Tool_type), "not_tool", reg.data.p1$Tool_type)
reg.data.p1$Tool_type = factor(reg.data.p1$Tool_type, levels = c(
  "not_tool", "notch", "denticulate", "notch/denticulate", "scraper", "point", "biface", "multiple", "other"
))

reg.data.p1$Weathering_class = factor(reg.data.p1$Weathering_class, levels = c("not_weathered", "strongly_weathered", "mildly_weathered", "weakly_weathered", "other"))
reg.data.p1$Artifact_type = factor(reg.data.p1$Artifact_type, levels = c("shatter", "complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment"))

fit.p1 = glm(recycled ~ ., family = binomial(), data = reg.data.p1)
summary(fit.p1) 
p1.df = tidy(fit.p1)
p1.df$location = "Semizbugu P1"
p1.df$signf = ifelse(p1.df$p.value < 0.05, "signf", "not signf")


reg.data.p2 = all_artifacts %>% filter(location == "Semizbugu P2") %>%
  dplyr::select("recycled", "Weathering_class", "Artifact_type", "Tool_type", "Cortex_percentage",
                "Length", "Width", "Thickness", 
                "Weight")
reg.data.p2$Tool_type = as.character(reg.data.p2$Tool_type)
reg.data.p2$Tool_type = ifelse(is.na(reg.data.p2$Tool_type), "not_tool", reg.data.p2$Tool_type)
reg.data.p2$Tool_type = factor(reg.data.p2$Tool_type, levels = c(
  "not_tool", "notch", "denticulate", "notch/denticulate", "scraper", "point", "biface", "multiple", "other"
))
reg.data.p2$Weathering_class = factor(reg.data.p2$Weathering_class, levels = c("not_weathered", "strongly_weathered", "mildly_weathered", "weakly_weathered", "other"))
reg.data.p2$Artifact_type = factor(reg.data.p2$Artifact_type, levels = c("shatter", "complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment"))

fit.p2 = glm(recycled ~ ., family = binomial(), data = reg.data.p2)
summary(fit.p2) 
p2.df = tidy(fit.p2)
p2.df$location = "Semizbugu P2"
p2.df$signf = ifelse(p2.df$p.value < 0.05, "signf", "not signf")


reg.data.p5 = all_artifacts %>% filter(location == "Semizbugu P5") %>%
  dplyr::select("recycled", "Weathering_class", "Artifact_type", "Tool_type", "Cortex_percentage",
                "Length", "Width", "Thickness", 
                "Weight")
reg.data.p5$Tool_type = as.character(reg.data.p5$Tool_type)
reg.data.p5$Tool_type = ifelse(is.na(reg.data.p5$Tool_type), "not_tool", reg.data.p5$Tool_type)
reg.data.p5$Tool_type = factor(reg.data.p5$Tool_type, levels = c(
  "not_tool", "notch", "denticulate", "notch/denticulate", "scraper", "point", "biface", "multiple", "other"
))
reg.data.p5$Weathering_class = factor(reg.data.p5$Weathering_class, levels = c("not_weathered", "strongly_weathered", "mildly_weathered", "weakly_weathered", "other"))
reg.data.p5$Artifact_type = factor(reg.data.p5$Artifact_type, levels = c("shatter", "complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment"))

fit.p5 = glm(recycled ~ ., family = binomial(), data = reg.data.p5)
summary(fit.p5) 
p5.df = tidy(fit.p5)
p5.df$location = "Semizbugu P5"
p5.df$signf = ifelse(p5.df$p.value < 0.05, "signf", "not signf")


reg.data.s10a = all_artifacts %>% filter(location == "Semizbugu 10A") %>%
  dplyr::select("recycled", "Weathering_class", "Artifact_type", "Tool_type", "Cortex_percentage",
                "Length", "Width", "Thickness", 
                "Weight")
reg.data.s10a$Tool_type = as.character(reg.data.s10a$Tool_type)
reg.data.s10a$Tool_type = ifelse(is.na(reg.data.s10a$Tool_type), "not_tool", reg.data.s10a$Tool_type)
reg.data.s10a$Tool_type = factor(reg.data.s10a$Tool_type, levels = c(
  "not_tool", "notch", "denticulate", "notch/denticulate", "scraper", "point", "biface", "multiple", "other"
))
reg.data.s10a$Weathering_class = factor(reg.data.s10a$Weathering_class, levels = c("not_weathered", "strongly_weathered", "mildly_weathered", "weakly_weathered", "other"))
reg.data.s10a$Artifact_type = factor(reg.data.s10a$Artifact_type, levels = c("shatter", "complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment"))

fit.s10a = glm(recycled ~ ., family = binomial(), data = reg.data.s10a)
summary(fit.s10a) 
s10a.df = tidy(fit.s10a)
s10a.df$location = "Semizbugu 10A"
s10a.df$signf = ifelse(s10a.df$p.value < 0.05, "signf", "not signf")


reg.data.s4 = all_artifacts %>% filter(location == "Semizbugu 4") %>%
  dplyr::select("recycled", "Weathering_class", "Artifact_type", "Tool_type", "Cortex_percentage",
                "Length", "Width", "Thickness", 
                "Weight")
reg.data.s4$Tool_type = as.character(reg.data.s4$Tool_type)
reg.data.s4$Tool_type = ifelse(is.na(reg.data.s4$Tool_type), "not_tool", reg.data.s4$Tool_type)
reg.data.s4$Tool_type = factor(reg.data.s4$Tool_type, levels = c(
  "not_tool", "notch", "denticulate", "notch/denticulate", "scraper", "point", "biface", "multiple", "other"
))
reg.data.s4$Weathering_class = factor(reg.data.s4$Weathering_class, levels = c("not_weathered", "strongly_weathered", "mildly_weathered", "weakly_weathered", "other"))
reg.data.s4$Artifact_type = factor(reg.data.s4$Artifact_type, levels = c("shatter", "complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment"))

fit.s4 = glm(recycled ~ ., family = binomial(), data = reg.data.s4)
summary(fit.s4) 
s4.df = tidy(fit.s4)
s4.df$location = "Semizbugu 4"
s4.df$signf = ifelse(s4.df$p.value < 0.05, "signf", "not signf")

#jtools::plot_summs(fit.p1, fit.p2, fit.p5, fit.s10a, fit.s4)

allregs = rbind(p1.df, p2.df, p5.df, s10a.df, s4.df)
allregs$lower = allregs$estimate - allregs$std.error
allregs$upper = allregs$estimate + allregs$std.error

allregs$term = factor(allregs$term, levels = 
                        c( "(Intercept)",
                           "Artifact_typecomplete_flake", "Artifact_typebroken_flake", "Artifact_typetool", "Artifact_typetool_fragment", "Artifact_typecore", "Artifact_typecore_fragment",
                           "Tool_typenotch", "Tool_typedenticulate", "Tool_typenotch/denticulate", "Tool_typescraper", "Tool_typepoint", "Tool_typebiface", "Tool_typemultiple", "Tool_typeother",
                           "Weathering_classstrongly_weathered", "Weathering_classmildly_weathered", "Weathering_classweakly_weathered", "Weathering_classother",
                           "Cortex_percentage", 
                           "Length", "Width", "Thickness", "Weight"
                        ))

regsp = ggplot(data = allregs %>% filter(signf == "signf")) +
  geom_hline(yintercept = 0, color = I("red")) +
  geom_linerange(aes(x=term, ymin = lower, ymax = upper, group = location, color = location), position = position_dodge2(width = 0.75)) +
  geom_point(aes(x=term, y = estimate, group = location, color = location), position = position_dodge2(width = 0.75), size = 4) +
  coord_flip() +
  scale_color_colorblind() +
  scale_x_discrete(labels = rev(c("weight", "thickness", "cortex percentage", "weathering class: other", "weathering class: weak", "weathering class: mild", "weathering class: strong", "tool type: multiple", "tool type: biface", "tool type: scraper", "artifact type: tool")))
plot(regsp)
ggsave(filename= "figures/regression-results-by-location.tiff", regsp, 
       dpi = 300, width = 8, height = 5)


