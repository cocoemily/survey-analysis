#Cortex Ratio calculation script for paleocore-sss collection forms
## units are millimeters (length, width, thickness) and grams (weight)

library(tidyverse)

artifacts = read_csv("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/cleaned_june_artifacts.csv")
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

calculate_cortex_ratio(cr.artifacts %>% filter(Site_name %in% c("Square 1", "Square 2", "Square 3")))
calculate_cortex_ratio(cr.artifacts %>% filter(Site_name %in% c("Square 4", "Square 5")))
calculate_cortex_ratio(cr.artifacts %>% filter(Site_name %in% c("Square 6")))
