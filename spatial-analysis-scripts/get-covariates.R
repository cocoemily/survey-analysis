library(tidyverse)
library(rgdal)
library(raster)
library(sf)
library(sfhotspot)
library(spDataLarge)
library(spdep)
library(tmap)
library(spatstat)

art.dens = density(ppp, sigma = bw.diggle)
plot(art.dens, useRaster = F, main = "Artifact density")
contour(art.dens, add = TRUE)

#### create covariate density plots ####
cf.ppp = as.ppp(st_data %>% filter(compl_flk == 1))
marks(cf.ppp) = NULL
Window(cf.ppp) = win
cf.dens = density(cf.ppp)
plot(cf.dens, useRaster = F)
contour(cf.dens, add = T)

bf.ppp = as.ppp(st_data %>% filter(broke_flk == 1))
marks(bf.ppp) = NULL
Window(bf.ppp) = win
bf.dens = density(bf.ppp)

tl.ppp = as.ppp(st_data %>% filter(tool == 1))
marks(tl.ppp) = NULL
Window(tl.ppp) = win
tl.dens = density(tl.ppp)

tlf.ppp = as.ppp(st_data %>% filter(tool_frag == 1))
marks(tlf.ppp) = NULL
Window(tlf.ppp) = win
tlf.dens = density(tlf.ppp)

cr.ppp = as.ppp(st_data %>% filter(core == 1))
marks(cr.ppp) = NULL
Window(cr.ppp) = win
cr.dens = density(cr.ppp)

crf.ppp = as.ppp(st_data %>% filter(core_frag == 1))
marks(crf.ppp) = NULL
Window(crf.ppp) = win
crf.dens = density(crf.ppp)

sh.ppp = as.ppp(st_data %>% filter(shatter == 1))
marks(sh.ppp) = NULL
Window(sh.ppp) = win
sh.dens = density(sh.ppp)
plot(sh.dens)

