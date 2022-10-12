library(tidyverse)
library(rgdal)
library(raster)
library(sf)
library(sfhotspot)
library(spDataLarge)
library(spdep)
library(tmap)
library(spatstat)

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


#### SPATIAL POINT PATTERN ANALYSIS IN R TESTING ####
st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))

ppp = as.ppp(st_data)
marks(ppp) = NULL
Window(ppp) = win

plot(ppp)
summary(ppp)

intensity(ppp) ##0.8 points per square unit

bw.diggle(ppp) #statistically best performance for kernel estimates
D = density(ppp, diggle = T,  sigma = bw.diggle, adjust = 2) #but choosing sigma with diggle assumes a Cox process (Poisson process with random intensity)
plot(D)

#### POINT PROCESS ANALYSIS ####
#### ppp with recycled artifacts ####
rcycl.ppp = as.ppp(st_data %>% filter(rcycl == 1))
marks(rcycl.ppp) = NULL
Window(rcycl.ppp) = win
plot(rcycl.ppp)

K2 = density(rcycl.ppp)
plot(K2, useRaster = F, main = "Recycled artifact density")
contour(K2, add = TRUE)

#test for complete spatial randomness
plot(Gest(rcycl.ppp))

plot(Fest(rcycl.ppp))
plot(Jest(rcycl.ppp))

Kest = Kest(rcycl.ppp)
plot(Kest)
#all measures confirm CSR --> homogenous process

#### recycling point process and artifact density ####
cdf.test(rcycl.ppp, art.dens, test = "ks")
#p value less than 0.05, so reject the null -> recycled points do depend on underlying density of artifacts







rho = rhohat(rcycl.ppp, art.dens)
plot(rho, las=1, main=NULL)
#increasing intensity of recycled artifacts as density of artifacts increases





#### MARKED POINT PROCESSES ####
atype.ppp = as.ppp(st_data %>% select(Artfct_t)) 
Window(atype.ppp) = win

summary(atype.ppp)
plot(split(atype.ppp))

#independ intensities
plot(density(split(atype.ppp)), useRaster = F)

#relative proportions of intensity ?? -- doesn't work right now
# Y = density(split(atype.ppp))
# attach(Y)
# pCF = eval.im(complete_flake/(complete_flake + broken_flake + tool + tool_fragment +
#                           + core + core_fragement + shatter))
# detach(Y)


plot(alltypes(atype.ppp, "Kdot"))

##need to determine which point process is best for modeling
## also need to figure out how to interpret the results 
test = ppm(atype.ppp, ~marks * K1)
summary(test)



