set.seed(120109)

library(tidyverse)
library(rgdal)
library(raster)
library(sf)
library(spatstat)

p1.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p1-artifacts")
p1.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p1-sqs")

data = p1.artifacts
window = p1.window

source("spatial-analysis-scripts/clean-data.R")

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))


artifact.ppp = as.ppp(st_data)
marks(artifact.ppp) = NULL
Window(artifact.ppp) = win
artifact.dens = density(artifact.ppp, sigma = bw.diggle)
#plot(artifact.dens)

rcycl.ppp = as.ppp(st_data %>% filter(rcycl == 1))
marks(rcycl.ppp) = NULL
Window(rcycl.ppp) = win
rcycl.dens3 = density(rcycl.ppp, sigma = bw.diggle)
plot(rcycl.dens3, useRaster = F)
plot(rcycl.ppp, add = T)

rcycl.ppm = ppm(rcycl.ppp ~ log(artifact.dens))
rel.rcycl.ppm = predict(rcycl.ppm)
plot(rel.rcycl.ppm, useRaster = F)
#plot(rcycl.ppm, add = T)
rcycl.rast = raster(rel.rcycl.ppm)
plot(rcycl.rast)

### testing with tools
tool.ppp = as.ppp(st_data %>% filter(tool == 1))
marks(tool.ppp) = NULL
Window(tool.ppp) = win
tool.dens = density(tool.ppp, sigma = bw.diggle)
plot(tool.dens3, useRaster = F)

tool.ppm = ppm(tool.ppp ~ log(artifact.dens))
rel.tool.ppm = predict(tool.ppm)
plot(rel.tool.ppm, useRaster = F)
tool.rast = raster(rel.tool.ppm)
plot(tool.rast)

plot(overlay(tool.rast, rcycl.rast, fun = function(r1, r2){return(r1-r2)}))
