library(tidyverse)
library(rgdal)
library(raster)
library(sf)
library(sfhotspot)
library(spDataLarge)
library(spdep)
library(tmap)
library(spatstat)

#### MARKED MULTITYPE POINT PROCESSES ####
atype.ppp = as.ppp(st_data %>% select(Artfct_t)) 
Window(atype.ppp) = win

summary(atype.ppp)
plot(split(atype.ppp))

#independent intensities
atype.dens = density(split(atype.ppp), sigma = bw.diggle, adjust = 2)
plot(atype.dens, useRaster = F)
plot(Smooth(atype.ppp), useRaster = F)



