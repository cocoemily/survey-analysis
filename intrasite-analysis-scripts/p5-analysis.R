## Poisson point process analysis for Semizbugu P5
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

p5.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p5-artifacts")
p5.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p5-sqs")

data = p5.artifacts
window = p5.window

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

Kin = Kinhom(ppp, lambda = Dnn, correction = "Ripley")
plot(Kin)
#Lin = Linhom(ppp, lambda = Dnn)
#plot(Lin)
#both K function falls slightly below the poisson line --> points are slightly more dispersed than expected
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
wdata = st_data %>% filter(Wthrng_c != "other" & !is.na(Wthrng_c))
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
segregation.test(ttppp, nsim=99)
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
##recycling ppp is inhomogeneous

#test for complete spatial randomness
mad.test(rcycl.ppp, Kinhom, nsims = 99, use.theo = T)
#mad test indicates no CSR
dclf.test(rcycl.ppp, Kinhom, nsims = 99, use.theo = T)
#dclf test indicate no CSR
hopskel.test(rcycl.ppp)
#hopskel test indicates no CSR

D = density(rcycl.ppp, sigma=bw.diggle)
plot(D, useRaster=F)

sqrt(nrow(rcycl.data))
Dnn = nndensity(rcycl.ppp, k = 24)
plot(Dnn, useRaster=F)

Ks = Kinhom(rcycl.ppp, lambda = Dnn)
plot(Ks)
Ks = Kinhom(rcycl.ppp, lambda = D)
plot(Ks)
Ls = Linhom(rcycl.ppp, lambda = Dnn)
plot(Ls)
Ls = Linhom(rcycl.ppp, lambda = D)
plot(Ls)
#both K and L fall slightly below the poisson line --> points are slightly more dispersed than expected
##model recycled PPP as Gibbs model, maybe a Poisson though?


##### Dependence of recycling points intensity on elevation and slope #####
##cdf null hypothesis - CDF of the covariate at all points is equal 
## to the CDF of covariate evaluated at the location of the point pattern

cdf.test(rcycl.ppp, dem_im)
cdf.test(rcycl.ppp, sp_im)
cdf.test(rcycl.ppp, sd_im)

auc(rcycl.ppp, dem_im)
auc(rcycl.ppp, sd_im)
auc(rcycl.ppp, sp_im)

##### Dependence of recycling points intensity on underlying artifact density ####
sqrt(nrow(st_data))
Dnn = nndensity(ppp, k = 40)
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

##### Gibbs process models -- TODO####


#### Marked point process for segregation analysis ####
rcycl.data = st_data %>% filter(rcycl == 1)
rcycl.data = rcycl.data %>% filter(poss_roll == F | is.na(poss_roll))
rcycl.ppp = as.ppp(rcycl.data)
marks(rcycl.ppp) = rcycl.data$Artfct_t
Window(rcycl.ppp) = win
plot(rcycl.ppp)
summary(rcycl.ppp)

plot(density(split(rcycl.ppp)), useRaster = F)


segregation.test(rcycl.ppp, nsim = 99) 
#no spatial variation in recycled object type

ttdata = rcycl.data %>% filter(!is.na(tool.type))
ttrpp = as.ppp(ttdata)
marks(ttrpp) = ttdata$tool.type
Window(ttrpp) = win
plot(ttrpp)
segregation.test(ttrpp, nsim = 99) ###not working??
plot(density(split(ttrpp)), useRaster = F)

#weathering class
wdata = rcycl.data %>% filter(!is.na(Wthrng_c))
wppp = as.ppp(wdata)
marks(wppp) = wdata$Wthrng_c
Window(wppp) = win
plot(wppp)
segregation.test(wppp, nsim=99)
#there is spatial segregation of weathering types among recycled implements
plot(density(split(wppp)), useRaster=F)

#### Figure -- recycling types spatial ####
sp.p5 = readOGR("data/artifact-shapefiles", "p5-artifacts")
sp.p5 = sp.p5[sp.p5$Id_nm %in% p5$Id_number,]

sp.p5@data = sp.p5@data %>%
  mutate(Recycling_type =
           ifelse(str_detect(Rcyclng_n, "none"), "none", 
                  ifelse(str_detect(Rcyclng_n, " "), "multiple", 
                         Rcyclng_n))) %>%
  mutate(Recycling_type = str_replace_all(Recycling_type, "_", " ")) %>%
  mutate(Recycling_type = ifelse(Recycling_type == "core on flake blade", "core on flake/blade", Recycling_type))
sp.p5@data$Recycling_type = factor(sp.p1@data$Recycling_type, 
                                   levels = c("double patina", "core on flake/blade", "core on hammerstone", "multiple", "other", "none"))

sf.p5 = st_as_sf(sp.p5)

pal1 = c("double patina" = "#CC6677", 
         "core on flake/blade" = "#DDCC77", 
         "core on hammerstone" = "#117733", 
         "multiple" = "#332288", 
         "other" = "#AA4499", 
         "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

plot.sf.p5 = sf.p5[!is.na(sf.p5$Recycling_type),]

rt.spat.plot = ggplot() +
  geom_sf(data = plot.sf.p5[plot.sf.p5$Recycling_type == "none",], size = 0.25, alpha = 0.25) +
  geom_sf(data = plot.sf.p5[plot.sf.p5$Recycling_type != "none",], aes(color = Recycling_type), size = 1) +
  #scale_color_colorblind() +
  scale_color_manual(values = pal1) +
  labs(color = "Recycling signature")
#ggsave(plot = rt.spat.plot, file = "~/Desktop/NYU/Dissertation-Research/papers/Coco_AK/figures/p1-recycling-signature-spatial.tiff")

plot(rt.spat.plot)

