set.seed(120109)

library(tidyverse)
library(rgdal)
library(raster)
library(sf)
library(spatstat)

p2.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p2-artifacts")
p2.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p2-sqs")

data = p2.artifacts
window = p2.window

source("spatial-analysis-scripts/clean-data.R")
#source("spatial-analysis-scripts/get-covariates.R")

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))


artifact.ppp = as.ppp(st_data)
marks(artifact.ppp) = NULL
Window(artifact.ppp) = win

rcycl.ppp = as.ppp(st_data %>% filter(rcycl == 1))
marks(rcycl.ppp) = NULL
Window(rcycl.ppp) = win


#### Hotspots Exploration ####
LR = scanLRTS(rcycl.ppp, r = 2*bw.diggle(artifact.ppp))
pvals = eval.im(pchisq(LR, df = 1, lower.tail=F))
rcycl.hs = as.im(pvals < 0.05)
Window(rcycl.hs) = win
plot(rcycl.hs)

aLR = scanLRTS(artifact.ppp, r = 2*bw.diggle(artifact.ppp))
aPvals = eval.im(pchisq(aLR, df = 1, lower.tail=F))
artifact.hs = as.im(aPvals < 0.05)
Window(artifact.hs) = win
plot(artifact.hs)

cfLR = scanLRTS(compl_flk.ppp, r = 2*bw.diggle(artifact.ppp))
cfPvals = eval.im(pchisq(cfLR, df = 1, lower.tail=F))
cf.hs = as.im(cfPvals < 0.05)
Window(cf.hs) = win
plot(cf.hs)

bfLR = scanLRTS(broke_flk.ppp, r = 2*bw.diggle(artifact.ppp))
bfPvals = eval.im(pchisq(bfLR, df = 1, lower.tail=F))
bf.hs = as.im(bfPvals < 0.05)
Window(bf.hs) = win
plot(bf.hs)

pppdata = st_data %>% filter(core == 1 | core_frag == 1)
allcore.ppp = as.ppp(pppdata)
marks(allcore.ppp) = NULL
Window(allcore.ppp) = win
coreLR = scanLRTS(allcore.ppp, r = 2*bw.diggle(artifact.ppp))
corePvals = eval.im(pchisq(coreLR, df = 1, lower.tail=F))
core.hs = as.im(corePvals < 0.05)
Window(core.hs) = win
plot(core.hs)

pppdata = st_data %>% filter(tool == 1 | tool_frag == 1)
alltool.ppp = as.ppp(pppdata)
marks(alltool.ppp) = NULL
Window(alltool.ppp) = win
toolLR = scanLRTS(alltool.ppp, r = 2*bw.diggle(artifact.ppp))
toolPvals = eval.im(pchisq(toolLR, df = 1, lower.tail=F))
tool.hs = as.im(toolPvals < 0.05)
Window(tool.hs) = win
plot(tool.hs)

plot(
  (artifact.hs*2) - rcycl.hs, 
  useRaster=F
)

plot(
  (cf.hs*2) - rcycl.hs, 
  useRaster=F
)

plot(
  (bf.hs*2) - rcycl.hs, 
  useRaster=F
)

plot(
  (core.hs*2) - rcycl.hs, 
  useRaster=F
)

plot(
  (tool.hs*2) - rcycl.hs, 
  useRaster=F
)
