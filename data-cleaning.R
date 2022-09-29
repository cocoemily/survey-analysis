#DATA CLEANING
library(tidyverse)
library(sp)
library(rgdal)

artifacts1 = read_csv("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/cleaned_june_artifacts.csv")
artifacts2 = read_csv("~/Desktop/NYU/Dissertation-Research/Survey/July-survey/cleaned_july_artifacts.csv")

paleocore_sss_artifact_form_all_versions_False_2022_08_01_05_36_18 <- read_delim("~/Desktop/NYU/Dissertation-Research/Survey/July-survey/paleocore-sss_artifact-form_-_all_versions_-_False_-_2022-08-01-05-36-18.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

s10 = paleocore_sss_artifact_form_all_versions_False_2022_08_01_05_36_18 %>% filter(Site_name == "Semizbugu 10A")
s4 = paleocore_sss_artifact_form_all_versions_False_2022_08_01_05_36_18 %>% filter(Site_name == "Semizbugu 4")

artifacts = rbind(artifacts1, artifacts2)
collections = rbind(s10, s4)

artifacts = artifacts %>% mutate(location = ifelse(Site_name %in% c("Square 1", "Square  1", "Square 2", "Square 3"), "Semizbugu P1", 
                                                   ifelse(Site_name %in% c("Square 4", "Square 5"), "Semizbugu P2", "Semizbugu P5")))
artifacts$recycled = !is.na(artifacts$Recycling_description)
artifacts$double_patina = str_detect(artifacts$Recycling_indications, "double_patina")
artifacts$Raw_material_description = tolower(artifacts$Raw_material_description)

artifacts = artifacts %>% filter(is.na(Problem_notes)) 
artifacts = artifacts %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers

collections$location = collections$Site_name
collections$recycled = !is.na(collections$Recycling_description)
collections$double_patina = str_detect(collections$Recycling_indications, "double_patina")
collections$Raw_material_description = tolower(collections$Raw_material_description)

collections = collections %>% filter(is.na(Problem_notes))
collections = collections %>%
  mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
         Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
         Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
  filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers

p1 = artifacts %>% filter(location == "Semizbugu P1")
p2 = artifacts %>% filter(location == "Semizbugu P2")
p5 = artifacts %>% filter(location == "Semizbugu P5")
s10a = collections %>% filter(location == 'Semizbugu 10A')
s4 = collections %>% filter(location == "Semizbugu 4")

#### WRITE CSV ####
write_csv(p1, file = "/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-dataframes/p1-artifacts.csv")
write_csv(p2, file = "/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-dataframes/p2-artifacts.csv")
write_csv(p5, file = "/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-dataframes/p5-artifacts.csv")
write_csv(s10a, file = "/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-dataframes/10a-artifacts.csv")
write_csv(s4, file = "/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-dataframes/4-artifacts.csv")

#### WRITE SHAPEFILES ####
p1_shp = SpatialPointsDataFrame(p1[,93:94], p1, proj4string = CRS("+proj=longlat +datum=WGS84"))
writeOGR(p1_shp, "/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles",
         "p1-artifacts", driver="ESRI Shapefile")

p2_shp = SpatialPointsDataFrame(p2[,93:94], p2, proj4string = CRS("+proj=longlat +datum=WGS84"))
writeOGR(p2_shp, "/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles",
         "p2-artifacts", driver="ESRI Shapefile")

p5_shp = SpatialPointsDataFrame(p5[,93:94], p5, proj4string = CRS("+proj=longlat +datum=WGS84"))
writeOGR(p5_shp, "/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles",
         "p5-artifacts", driver="ESRI Shapefile")
