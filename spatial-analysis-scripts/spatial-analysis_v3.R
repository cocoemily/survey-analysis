#non-Poisson models

library(tidyverse)
library(rgdal)
library(raster)
library(sf)
library(sfhotspot)
library(spDataLarge)
library(spdep)
library(tmap)
library(spatstat)
library(AICcmodavg)

set.seed(120109)

p1.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p1-artifacts")
p2.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p2-artifacts")
p5.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p5-artifacts")

p1.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/p1-sqs.shp", layer = "p1-sqs")
p2.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/p1-sqs.shp", layer = "p2-sqs")
p5.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/p1-sqs.shp", layer = "p5-sqs")


data = p1.artifacts
window = p1.window

source("spatial-analysis-scripts/clean-data.R")
source("spatial-analysis-scripts/get-covariates.R")

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))

rcycl.ppp = as.ppp(st_data %>% filter(rcycl == 1))
marks(rcycl.ppp) = NULL
Window(rcycl.ppp) = win
plot(rcycl.ppp)


#Cox models -recycled point pattern is randomly influenced by unobserved covariates
kppm1 = kppm(rcycl.ppp ~ 1, "Thomas")
kppm1
plot(kppm1)


#Gibbs model - interactions occur between points of the recycled point pattern
gppm1 = ppm(rcycl.ppp ~ 1, Strauss(5))
gppm1

