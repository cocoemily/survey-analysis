#Cortex Ratio calculation script for paleocore-sss collection forms
## units are millimeters (length, width, thickness) and grams (weight)

library(tidyverse)

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


all_artifacts = rbind(
  p1 %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight, Cortex_percentage, Cortex_description), 
  p2 %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight,Cortex_percentage, Cortex_description), 
  p5 %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight, Cortex_percentage, Cortex_description), 
  s10a %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight, Cortex_percentage, Cortex_description), 
  s4 %>% select(Id_number, location, recycled, Raw_material_description, Weathering_class, Artifact_type, Flake_thickness, Flake_length, Flake_width, Maximum_core_length, Maximum_core_width, Maximum_core_thickness, Weight, Cortex_percentage, Cortex_description)
)
cr.artifacts = artifacts %>% filter(!is.na(Cortex_percentage))

raw.mat = unique(tolower(cr.artifacts$Raw_material_description))

#calculations without core thickness measurement
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
  bifaces$surface.area = (
    4*pi*(((((bifaces$Flake_length/2) ^ 1.6075)*((bifaces$Flake_width/2) ^1.6075)+
              ((bifaces$Flake_length/2) ^1.6075)*((bifaces$Flake_thickness/2) ^1.6075)+
              ((bifaces$Flake_width/2)^1.6075)*((bifaces$Flake_thickness/2)^1.6075))/3)^(1/1.6075))
  )
  bifaces$cortex.area = bifaces$surface.area*(bifaces$Cortex_percentage/100.0)
  bifaces = bifaces %>% filter(!is.na(cortex.area))
  
  ## Core surface area and cortex area calculation (mm2)
  #cores$alpha = acos((pmin(cores$Maximum_core_length, cores$Maximum_core_width)/pmax(cores$Maximum_core_length, cores$Maximum_core_width)))
  cores$surface.area = (
    4*pi*(((((cores$Maximum_core_length/2) ^ 1.6075)*((cores$Maximum_core_width/2) ^1.6075)+
              ((cores$Maximum_core_length/2) ^1.6075)*((cores$Maximum_core_thickness/2) ^1.6075)+
              ((cores$Maximum_core_width/2)^1.6075)*((cores$Maximum_core_thickness/2)^1.6075))/3)^(1/1.6075))
  )
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
    summarize(count = n())
  nodule.count = sum(core.count$count)#need to think about how this is being calculated
  nodule.volume = assemblage.volume/nodule.count
  
  ## Theoretical surface area calculations
  theor.sa = (4*pi*((3*nodule.volume)/(4*pi))^(2/3)) #check to see what model for nodule volume 
  te.cortical.sa = theor.sa*nodule.count
  
  ## Cortex ratio calculation
  cr = to.cortical.sa / te.cortical.sa
  
  return(cr)
  
}

calculate_cortex_ratio(cr.artifacts %>% filter(location == "Semizbugu P1"))
calculate_cortex_ratio(cr.artifacts %>% filter(location == "Semizbugu P2"))
calculate_cortex_ratio(cr.artifacts %>% filter(location == "Semizbugu P5"))
#calculate_cortex_ratio(cr.artifacts %>% filter(location == "Semizbugu 10A"))
#calculate_cortex_ratio(cr.artifacts %>% filter(location == "Semizbugu 4"))
