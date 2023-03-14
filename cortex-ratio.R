## units are millimeters (length, width, thickness) and grams (weight)

library(tidyverse)

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


####CORTEX CLEANING####
cr.artifacts = all_artifacts %>% filter(!is.na(Cortex_percentage))
raw.mat = unique(tolower(cr.artifacts$Raw_material_description))


####Cortex Ratio calculations####
calculate_cortex_ratio = function(data) {
  flakes = data %>% filter(str_detect(Artifact_type, "flake"))
  cores = data %>% filter(str_detect(Artifact_type, "core"))
  tools = data %>% filter(str_detect(Artifact_type, "tool")) %>% filter(Tool_type != "biface")
  bifaces = data %>% filter(str_detect(Artifact_type, "tool")) %>% filter(Tool_type == "biface")
  
  
  ## Flake and Tool surface area and cortex area calcuation (mm2)
  flakes$surface.area = flakes$Flake_length*flakes$Flake_width
  flakes$cortex.area = flakes$surface.area*(flakes$Cortex_percentage/100.0)
  flakes = flakes %>% filter(!is.na(cortex.area))
  
  tools$surface.area = tools$Flake_length*tools$Flake_width
  tools$cortex.area = tools$surface.area*(tools$Cortex_percentage/100.0)
  tools = tools %>% filter(!is.na(cortex.area))
  
  ## Biface surface area and cortex area calculation (mm2)
  #surface area of scalene ellipsoid
  bifaces$surface.area = (
    4*pi*(((((bifaces$Flake_length/2) ^ 1.6075)*((bifaces$Flake_width/2) ^1.6075)+
              ((bifaces$Flake_length/2) ^1.6075)*((bifaces$Flake_thickness/2) ^1.6075)+
              ((bifaces$Flake_width/2)^1.6075)*((bifaces$Flake_thickness/2)^1.6075))/3)^(1/1.6075))
  )
  bifaces$cortex.area = bifaces$surface.area*(bifaces$Cortex_percentage/100.0)
  bifaces = bifaces %>% filter(!is.na(cortex.area))
  
  ## Core surface area and cortex area calculation (mm2)
  #surface area of scalene ellipsoid
  if(all(!is.na(cores$Maximum_core_thickness)) == TRUE){
    cores$surface.area = (
      4*pi*(((((cores$Maximum_core_length/2) ^ 1.6075)*((cores$Maximum_core_width/2) ^1.6075)+
                ((cores$Maximum_core_length/2) ^1.6075)*((cores$Maximum_core_thickness/2) ^1.6075)+
                ((cores$Maximum_core_width/2)^1.6075)*((cores$Maximum_core_thickness/2)^1.6075))/3)^(1/1.6075))
    )
  } else {
    cores$max = pmax(cores$Maximum_core_length, cores$Maximum_core_width)
    cores$surface.area = (
      4*pi*(((((cores$Maximum_core_length/2) ^ 1.6075)*((cores$Maximum_core_width/2) ^1.6075)+
                ((cores$Maximum_core_length/2) ^1.6075)*((cores$max/2) ^1.6075)+
                ((cores$Maximum_core_width/2)^1.6075)*((cores$max/2)^1.6075))/3)^(1/1.6075))
    )
  }
  
  cores$cortex.area = cores$surface.area*(cores$Cortex_percentage/100.0)
  cores = cores %>% filter(!is.na(cortex.area))
  
  ## Total observed cortical surface area calculation (mm2)
  to.cortical.sa = sum(flakes$cortex.area, na.rm = T) + sum(cores$cortex.area, na.rm = T) +
    sum(tools$cortex.area, na.rm = T) + sum(bifaces$cortex.area, na.rm = T)
  
  
  ## Assemblage weight calculation (g)
  assemblage.weight = sum(flakes$Weight, na.rm = T) + sum(cores$Weight, na.rm = T) +
    sum(tools$Weight, na.rm = T) + sum(bifaces$Weight, na.rm = T) 
  
  ## Assemblage volume calculation
  #look up rock density as necessary
  assemblage.volume <- assemblage.weight/0.00265 #weight divided by rock density
  
  ## Nodule volume calculation
  core.count  = cores %>% group_by(Raw_material_description) %>% 
    summarize(count = n(), 
              average.length = mean(Maximum_core_length), 
              average.width = mean(Maximum_core_width), 
              average.thick = mean(Maximum_core_thickness), 
              vol = average.length * average.width * average.thick, 
              total.vol = count * vol)
  est.cob.vol = sum(core.count$total.vol)
  #nodule.count = assemblage.volume/est.cob.vol
  nodule.count = sum(core.count$count)
  nodule.volume = assemblage.volume/nodule.count
  
  
  ## Theoretical surface area calculations
  #theor.sa = (4*pi*((3*nodule.volume)/(4*pi))^(2/3)) # check to see what model for nodule volume 
  theor.sa = 6*(nodule.volume^(2/3))
  te.cortical.sa = theor.sa*nodule.count
  
  ## Cortex ratio calculation
  cr = to.cortical.sa / te.cortical.sa
  
  return(cr)
  
}

calculate_cortex_ratio(all_artifacts %>% filter(location == "Semizbugu P1"))
calculate_cortex_ratio(all_artifacts %>% filter(location == "Semizbugu P2"))
calculate_cortex_ratio(all_artifacts %>% filter(location == "Semizbugu P5"))
calculate_cortex_ratio(all_artifacts %>% filter(location == "Semizbugu 10A"))
calculate_cortex_ratio(all_artifacts %>% filter(location == "Semizbugu 4"))

var(x = c(
  calculate_cortex_ratio(all_artifacts %>% filter(location == "Semizbugu P1")),
  calculate_cortex_ratio(all_artifacts %>% filter(location == "Semizbugu P2")), 
  calculate_cortex_ratio(all_artifacts %>% filter(location == "Semizbugu P5")), 
  calculate_cortex_ratio(all_artifacts %>% filter(location == "Semizbugu 4"))
))
