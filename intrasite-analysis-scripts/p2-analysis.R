## Poisson point process analysis for Semizbugu P2
set.seed(120109)

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

p2.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p2-artifacts")
p2.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p2-sqs")

data = p2.artifacts
window = p2.window

source("spatial-analysis-scripts/clean-data.R")
#source("spatial-analysis-scripts/get-covariates.R")

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))

####ALL POINTS####
quadrat.test(ppp, 5, method = "MonteCarlo") 
#results of quadrat test indicate a inhomogeneous process
plot(envelope(ppp, fun = Gest, nsim = 99))
plot(envelope(ppp, fun = Fest, nsim = 99))
#confirms inhomogeneous process

plot(ppp)

D = density(ppp, sigma=bw.diggle)
plot(D)

sqrt(nrow(st_data))
Dnn = nndensity(ppp, k = 25)
plot(Dnn)

mad.test(ppp, Linhom, nsims = 99, use.theo = T)
dclf.test(ppp, Linhom, nsims = 99, use.theo = T)
#both tests show spatial dependence of points

Kin = Kinhom(ppp, lambda = Dnn)
plot(Kin)
Lin = Linhom(ppp, lambda = Dnn)
plot(Lin)
#both K and L fall slight below the poisson line --> points are slightly more dispersed than expected
## total point process needs to be assessed as a Gibbs model?
plot(dem_im)

cdf.test(ppp, dem_im)
cdf.test(ppp, sp_im)
cdf.test(ppp, sd_im)
#location of points dependent on elevation and slope

###### TODO####
fit1 = ppm(ppp ~ dem_im)
summary(fit1)
plot(fit1, useRaster=F)
fitted.artifact.dens = intensity(fit1)

###### segregation tests #####
#null = spatially constant
#artifact types
ppp = as.ppp(st_data)
marks(ppp) = st_data$Artfct_t
Window(ppp) = win
plot(ppp)
segregation.test(ppp, nsim=99)
#there is spatial segregation of artifact types
plot(density(split(ppp)), useRaster=F)

#weathering class with other weathering
ppp = as.ppp(st_data)
marks(ppp) = st_data$Wthrng_c
Window(ppp) = win
plot(ppp)
segregation.test(ppp, nsim=99) 
#there is spatial segregation of weathering types
plot(density(split(ppp)), useRaster=F)

#weathering class without other weathering
wdata = st_data %>% filter(Wthrng_c != "other")
wppp = as.ppp(wdata)
marks(wppp) = wdata$Wthrng_c
Window(wppp) = win
plot(wppp)
segregation.test(wppp, nsim=99)
plot(density(split(wppp)), useRaster=F)
#spatial segregation of weathering types

ttdata = st_data %>% filter(!is.na(tool.type))
ttppp = as.ppp(ttdata)
marks(ttppp) = ttdata$tool.type
Window(ttppp) = win
plot(ttppp)
segregation.test(ttppp, nsim=99) #not working??
plot(density(split(ttppp)), useRaster=F)


#### ppp with recycled artifacts ####
rcycl.data = st_data %>% filter(rcycl == 1)
rcycl.data = rcycl.data %>% filter(poss_roll == F | is.na(poss_roll))
rcycl.ppp = as.ppp(rcycl.data)
marks(rcycl.ppp) = NULL
Window(rcycl.ppp) = win
plot(rcycl.ppp)

#test for homogeneity
quadrat.test(rcycl.ppp, 5, method = "MonteCarlo") 
plot(envelope(rcycl.ppp, fun = Gest, nsim = 99))
plot(envelope(rcycl.ppp, fun = Fest, nsim = 99))
##recycling ppp is mostly homogeneous, except at larger distances

#test for complete spatial randomness
mad.test(rcycl.ppp, Kest, nsims = 99, use.theo = T)
#mad test indicates no CSR
dclf.test(rcycl.ppp, Kest, nsims = 99, use.theo = T)
#dclf test indicate no CSR
hopskel.test(rcycl.ppp)
#hopskel test indicates no CSR

D = density(rcycl.ppp, sigma=bw.diggle)
plot(D, useRaster=F)

sqrt(nrow(rcycl.data))
Dnn = nndensity(rcycl.ppp, k = 11)
plot(Dnn, useRaster=F)

Ks = Kest(rcycl.ppp, lambda = Dnn)
plot(Ks)
Ks = Kest(rcycl.ppp, lambda = D)
plot(Ks)
Ls = Lest(rcycl.ppp, lambda = Dnn)
plot(Ls)
Ls = Lest(rcycl.ppp, lambda = D)
plot(Ls)
#both K and L fall above the poisson line --> points are more clustered than expected
##model recycled PPP as Cox/Cluster models, maybe a Poisson though?


##### Dependence of recycling points intensity on elevation and slope #####
##cdf null hypothesis - CDF of the covariate at all points is equal 
## to the CDF of covariate evaluated at the location of the point pattern

cdf.test(rcycl.ppp, dem_im)
cdf.test(rcycl.ppp, sp_im)
cdf.test(rcycl.ppp, sd_im)

auc(rcycl.ppp, sd_im)
auc(rcycl.ppp, sp_im)

##### Dependence of recycling points intensity on underlying artifact density ####
sqrt(nrow(st_data))
Dnn = nndensity(ppp, k = 25)
plot(Dnn)

cdf.test(rcycl.ppp, as.im(Dnn))
berman.test(rcycl.ppp, as.im(Dnn))
berman.test(rcycl.ppp, as.im(Dnn), "Z2")

auc(rcycl.ppp, as.im(Dnn))


##### Dependence of recycling points intensity on underlying density of retouched artifacts ####
cdf.test(rcycl.ppp, as.im(retouch.dens))
berman.test(rcycl.ppp, as.im(retouch.dens))
berman.test(rcycl.ppp, as.im(retouch.dens), "Z2")

auc(rcycl.ppp, as.im(retouch.dens))

##### Volumetric and weight dependence #####
cdf.test(rcycl.ppp, as.im(weight.dens))
auc(rcycl.ppp, as.im(weight.dens))

cdf.test(rcycl.ppp, as.im(thick.dens))
auc(rcycl.ppp, as.im(thick.dens))

cdf.test(rcycl.ppp, as.im(length.dens))
auc(rcycl.ppp, as.im(length.dens))

cdf.test(rcycl.ppp, as.im(width.dens))
auc(rcycl.ppp, as.im(width.dens))

##### Cox process models -- NOT WORKING ####
fit0 = kppm(rcycl.ppp ~ 1, "LGCP", method = "clik2")
AIC(fit0)
fit1 = kppm(rcycl.ppp ~ artifact.dens, clusters = "LGCP", method = "clik2")
AIC(fit1)
summary(fit1)
plot(fit1, useRaster=F)

fit2 = kppm(rcycl.ppp ~ artifact.dens + retouch.dens + cortex.dens +
              compl_flk.dens + broke_flk.dens + tool.dens + tool_frag.dens +
              core.dens , clusters = "LGCP", method = "clik2")
AIC(fit2)
summary(fit2)

fit3 = kppm(rcycl.ppp ~ artifact.dens + 
              weight.dens + length.dens + width.dens + thick.dens +
              mid_weather.dens + weak_weather.dens + str_weather.dens +
              retouch.dens + cortex.dens +
              compl_flk.dens + broke_flk.dens + tool.dens + tool_frag.dens +
              core.dens + core_frag.dens
            , clusters = "LGCP", method = "clik2")
AIC(fit3)
drop1(fit3)

fit4 = kppm(rcycl.ppp ~ artifact.dens + 
              weight.dens + length.dens + width.dens + thick.dens
            , clusters = "LGCP", method = "clik2")
AIC(fit4)
drop1(fit4)
summary(fit4)

fit5 = kppm(rcycl.ppp ~ artifact.dens + 
              mid_weather.dens + weak_weather.dens + str_weather.dens 
            , clusters = "LGCP", method = "clik2")
AIC(fit5) #poor fit



#### Marked point process for segregation analysis ####
marks(rcycl.ppp) = rcycl.data$Artfct_t
plot(rcycl.ppp)
plot(unmark(rcycl.ppp))
summary(rcycl.ppp)

plot(density(split(rcycl.ppp)), useRaster = F)

lambda = intensity(rcycl.ppp)
probs = lambda/sum(lambda) #for homogeneous process

segregation.test(rcycl.ppp, nsim = 99) 
#no spatial variation in recycled object type

ttdata = rcycl.data %>% filter(!is.na(tool.type))
ttrpp = as.ppp(ttdata)
marks(ttrpp) = ttdata$tool.type
plot(ttrpp)
segregation.test(ttrpp, nsim = 99) 
#no segregration by tool type

marks(rcycl.ppp) = rcycl.data$Wthrng_c
plot(rcycl.ppp)
segregation.test(rcycl.ppp, nsim = 99) 
#there is spatial segregation in weathering of recycled objects
plot(density(split(rcycl.ppp)), useRaster = F)

#weathering class without other weathering
wdata = rcycl.data %>% filter(Wthrng_c != "other")
wppp = as.ppp(wdata)
marks(wppp) = wdata$Wthrng_c
Window(wppp) = win
plot(wppp)
segregation.test(wppp, nsim=99)
plot(density(split(wppp)), useRaster=F)
