library(tidyverse)
library(rgdal)
library(raster)
library(sf)
library(sfhotspot)
library(spDataLarge)
library(spdep)
library(tmap)
library(spatstat)

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))

ppp = as.ppp(st_data)
marks(ppp) = NULL
Window(ppp) = win
artifact.dens = density(ppp, sigma = bw.diggle, adjust = 2)
plot(artifact.dens, useRaster = F, main = "Artifact density")
contour(artifact.dens, add = TRUE)

#### create covariate real density plots with Diggle correction ####
compl_flk.ppp = as.ppp(st_data %>% filter(compl_flk == 1))
marks(compl_flk.ppp) = NULL
Window(compl_flk.ppp) = win
compl_flk.dens = density(compl_flk.ppp, sigma = bw.diggle, adjust = 2)
plot(compl_flk.dens, useRaster = F)
contour(compl_flk.dens, add = T)

broke_flk.ppp = as.ppp(st_data %>% filter(broke_flk == 1))
marks(broke_flk.ppp) = NULL
Window(broke_flk.ppp) = win
broke_flk.dens = density(broke_flk.ppp, sigma = bw.diggle, adjust = 2)

tool.ppp = as.ppp(st_data %>% filter(tool == 1))
marks(tool.ppp) = NULL
Window(tool.ppp) = win
tool.dens = density(tool.ppp, sigma = bw.diggle)

tool_frag.ppp = as.ppp(st_data %>% filter(tool_frag == 1))
marks(tool_frag.ppp) = NULL
Window(tool_frag.ppp) = win
tool_frag.dens = density(tool_frag.ppp, sigma = bw.diggle, adjust = 2)

core.ppp = as.ppp(st_data %>% filter(core == 1))
marks(core.ppp) = NULL
Window(core.ppp) = win
core.dens = density(core.ppp, sigma = bw.diggle, adjust = 2)

core_frag.ppp = as.ppp(st_data %>% filter(core_frag == 1))
marks(core_frag.ppp) = NULL
Window(core_frag.ppp) = win
core_frag.dens = density(core_frag.ppp, sigma = bw.diggle, adjust = 2)
plot(core_frag.dens)
plot(core_frag.ppp, add = T)

shatter.ppp = as.ppp(st_data %>% filter(shatter == 1))
marks(shatter.ppp) = NULL
Window(shatter.ppp) = win
shatter.dens = density(shatter.ppp, sigma = bw.diggle, adjust = 2)

str_weather.ppp = as.ppp(st_data %>% filter(str_weather == 1))
marks(str_weather.ppp) = NULL
Window(str_weather.ppp) = win
str_weather.dens = density(str_weather.ppp, sigma = bw.diggle, adjust = 2)

mid_weather.ppp = as.ppp(st_data %>% filter(mid_weather == 1))
marks(mid_weather.ppp) = NULL
Window(mid_weather.ppp) = win
mid_weather.dens = density(mid_weather.ppp, sigma = bw.diggle, adjust = 2)

weak_weather.ppp = as.ppp(st_data %>% filter(weak_weather == 1))
marks(weak_weather.ppp) = NULL
Window(weak_weather.ppp) = win
weak_weather.dens = density(weak_weather.ppp, sigma = bw.diggle, adjust = 2)

not_weather.ppp = as.ppp(st_data %>% filter(not_weather == 1))
marks(not_weather.ppp) = NULL
Window(not_weather.ppp) = win
not_weather.dens = density(not_weather.ppp, sigma = bw.diggle, adjust = 2)

oth_weather.ppp = as.ppp(st_data %>% filter(oth_weather == 1))
marks(oth_weather.ppp) = NULL
Window(oth_weather.ppp) = win
oth_weather.dens = density(oth_weather.ppp, sigma = bw.diggle, adjust = 2)

type_flake.ppp = as.ppp(st_data %>% filter(type_flake == 1))
marks(type_flake.ppp) = NULL
Window(type_flake.ppp) = win
type_flake.dens = density(type_flake.ppp, sigma = bw.diggle, adjust = 2)

type_blade.ppp = as.ppp(st_data %>% filter(type_blade == 1))
marks(type_blade.ppp) = NULL
Window(type_blade.ppp) = win
type_blade.dens = density(type_blade.ppp, sigma = bw.diggle, adjust = 2)

type_bladelet.ppp = as.ppp(st_data %>% filter(type_bladelet == 1))
marks(type_bladelet.ppp) = NULL
Window(type_bladelet.ppp) = win
type_bladelet.dens = density(type_bladelet.ppp, sigma = bw.diggle, adjust = 2)

type_oth.ppp = as.ppp(st_data %>% filter(type_oth == 1))
marks(type_oth.ppp) = NULL
Window(type_oth.ppp) = win
type_oth.dens = density(type_oth.ppp, sigma = bw.diggle, adjust = 2)

rside_dorsal.ppp = as.ppp(st_data %>% filter(rside_dorsal == 1))
marks(rside_dorsal.ppp) = NULL
Window(rside_dorsal.ppp) = win
rside_dorsal.dens = density(rside_dorsal.ppp, sigma = bw.diggle, adjust = 2)

rside_ventral.ppp = as.ppp(st_data %>% filter(rside_ventral == 1))
marks(rside_ventral.ppp) = NULL
Window(rside_ventral.ppp) = win
rside_ventral.dens = density(rside_ventral.ppp, sigma = bw.diggle, adjust = 2)

rside_bifacial.ppp = as.ppp(st_data %>% filter(rside_bifacial == 1))
marks(rside_bifacial.ppp) = NULL
Window(rside_bifacial.ppp) = win
rside_bifacial.dens = density(rside_bifacial.ppp, sigma = bw.diggle, adjust = 2)

ttype_notch.ppp = as.ppp(st_data %>% filter(ttype_notch == 1))
marks(ttype_notch.ppp) = NULL
Window(ttype_notch.ppp) = win
ttype_notch.dens = density(ttype_notch.ppp, sigma = bw.diggle, adjust = 2)

ttype_dent.ppp = as.ppp(st_data %>% filter(ttype_dent == 1))
marks(ttype_dent.ppp) = NULL
Window(ttype_dent.ppp) = win
ttype_dent.dens = density(ttype_dent.ppp, sigma = bw.diggle, adjust = 2)

ttype_nd.ppp = as.ppp(st_data %>% filter(ttype_nd == 1))
marks(ttype_nd.ppp) = NULL
Window(ttype_nd.ppp) = win
ttype_nd.dens = density(ttype_nd.ppp, sigma = bw.diggle, adjust = 2)

ttype_scraper.ppp = as.ppp(st_data %>% filter(ttype_scraper == 1))
marks(ttype_scraper.ppp) = NULL
Window(ttype_scraper.ppp) = win
ttype_scraper.dens = density(ttype_scraper.ppp, sigma = bw.diggle, adjust = 2)

ttype_point.ppp = as.ppp(st_data %>% filter(ttype_point == 1))
marks(ttype_point.ppp) = NULL
Window(ttype_point.ppp) = win
ttype_point.dens = density(ttype_point.ppp, sigma = bw.diggle, adjust = 2)

ttype_biface.ppp = as.ppp(st_data %>% filter(ttype_biface == 1))
marks(ttype_biface.ppp) = NULL
Window(ttype_biface.ppp) = win
ttype_biface.dens = density(ttype_biface.ppp, sigma = bw.diggle, adjust = 2)

ttype_mult.ppp = as.ppp(st_data %>% filter(ttype_mult == 1))
marks(ttype_mult.ppp) = NULL
Window(ttype_mult.ppp) = win
ttype_mult.dens = density(ttype_mult.ppp, sigma = bw.diggle, adjust = 2)

ttype_oth.ppp = as.ppp(st_data %>% filter(ttype_oth == 1))
marks(ttype_oth.ppp) = NULL
Window(ttype_oth.ppp) = win
ttype_oth.dens = density(ttype_oth.ppp, sigma = bw.diggle, adjust = 2)
