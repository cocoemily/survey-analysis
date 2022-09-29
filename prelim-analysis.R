library(tidyverse)

artifacts = read_csv("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/cleaned_june_artifacts.csv")

summary(artifacts$`RMS E`)
summary(artifacts$`RMS N`)
summary(artifacts$`Elev. RMS`)

##ARTIFACT TYPE
table(artifacts$Artifact_type)

table(artifacts$Tool_type)
table(artifacts$Core_type)
table(artifacts$Flake_type)

##ARTIFACT SIZE
summary(artifacts$Flake_length)
summary(artifacts$Flake_width)
summary(artifacts$Flake_thickness)

summary(artifacts$Maximum_core_length)
summary(artifacts$Maximum_core_width)
summary(artifacts$Maximum_core_thickness)


##RECYCLING
table(artifacts$Recycling_indications)

artifacts[which(str_detect(artifacts$Recycling_indications, "none")), ]$Recycling_indications = "none"
artifacts$recycled = !is.na(artifacts$Recycling_description)
artifacts$double_patina = str_detect(artifacts$Recycling_indications, "double_patina")

table(artifacts$double_patina)
table(artifacts$recycled)

##WEATHERING
table(artifacts$Weathering_class)


##Density
dens = artifacts %>% dplyr::group_by(Site_name) %>%
  summarize(total = n(), 
            density =
              n()/ifelse(Site_name %in% c("Square 1", "Square 4", "Square 5", "Square 6"), 
                                       (20*20), (10*20)))
            




