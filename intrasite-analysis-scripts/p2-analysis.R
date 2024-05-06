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
library(ggpubr)

theme_set(theme_bw())

p2.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p2-artifacts")
p2.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p2-sqs")

data = p2.artifacts
window = p2.window

source("intrasite-analysis-scripts/clean-data.R")
#source("spatial-analysis-scripts/get-covariates.R")

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))

####ALL POINTS####
# 
# 
# quadrat.test(ppp, 5, method = "MonteCarlo") 
# #results of quadrat test indicate a inhomogeneous process
# # plot(envelope(ppp, fun = Gest, nsim = 99))
# # plot(envelope(ppp, fun = Fest, nsim = 99))
# #confirms inhomogeneous process
# 
# plot(ppp)
# 
# Lfun = envelope(ppp, fun = Linhom, nsim = 99, verbose = F)
# plot(Lfun)
# Lfun.global = envelope(ppp, Linhom, nsim = 19, rank = 1, global = T)
# plot(Lfun.global)
# 
# D = density(ppp, sigma=bw.diggle)
# plot(D)
# 
# sqrt(nrow(st_data))
# Dnn = nndensity(ppp, k = 25)
# plot(Dnn)
# 
# mad.test(ppp, Linhom, nsims = 99, use.theo = T)
# #dclf.test(ppp, Linhom, nsims = 99, use.theo = T)
# #both tests show spatial dependence of points
# 
# Kin = Kinhom(ppp, lambda = Dnn)
# plot(Kin)
# Lin = Linhom(ppp, lambda = Dnn)
# plot(Lin)
# #both K and L fall slight below the poisson line --> points are slightly more dispersed than expected
# 
# 
# ###### segregation tests #####
# #null = spatially constant
# #artifact types
# ppp = as.ppp(st_data)
# marks(ppp) = factor(st_data$Artfct_t)
# Window(ppp) = win
# plot(ppp)
# segregation.test(ppp, nsim=99)
# #there is spatial segregation of artifact types
# plot(density(split(ppp)), useRaster=F)
# 
# #weathering class with other weathering
# ppp = as.ppp(st_data)
# marks(ppp) = factor(st_data$Wthrng_c)
# Window(ppp) = win
# plot(ppp)
# segregation.test(ppp, nsim=99) 
# #there is spatial segregation of weathering types
# plot(density(split(ppp)), useRaster=F)
# 
# #weathering class without other weathering
# wdata = st_data %>% filter(Wthrng_c != "other")
# wppp = as.ppp(wdata)
# marks(wppp) = wdata$Wthrng_c
# Window(wppp) = win
# plot(wppp)
# segregation.test(wppp, nsim=99)
# plot(density(split(wppp)), useRaster=F)
# #spatial segregation of weathering types
# 
# ttdata = st_data %>% filter(!is.na(tool.type)) %>% filter(tool.type != "notch/denticulate")
# ttdata$tool.type = droplevels(ttdata$tool.type)
# ttppp = as.ppp(ttdata)
# marks(ttppp) = as.factor(ttdata$tool.type)
# Window(ttppp) = win
# plot(ttppp)
# segregation.test(ttppp, nsim=99) #not working??
# plot(density(split(ttppp)), useRaster=F)


#### RECYCLED ARTIFACTS ####
rcycl.data = st_data %>% filter(rcycl == 1)
rcycl.data = rcycl.data %>% filter(poss_roll == F | is.na(poss_roll))
rcycl.ppp = as.ppp(rcycl.data)
marks(rcycl.ppp) = NULL
Window(rcycl.ppp) = win
plot(rcycl.ppp)

#test for homogeneity
quadrat.test(rcycl.ppp, 5, method = "MonteCarlo") 
# plot(envelope(rcycl.ppp, fun = Gest, nsim = 99))
# plot(envelope(rcycl.ppp, fun = Fest, nsim = 99))
##recycling ppp is mostly homogeneous, except at larger distances

#test for complete spatial randomness
mad.test(rcycl.ppp, Lest, nsims = 99, use.theo = T)
#mad test indicates no CSR
dclf.test(rcycl.ppp, Kest, nsims = 99, use.theo = T)
#dclf test indicate no CSR
hopskel.test(rcycl.ppp)
#hopskel test indicates no CSR

rLfun = envelope(rcycl.ppp, fun = Lest, nsim = 99, verbose = F, correction = "Ripley")
plot(rLfun)
rLfun.global = envelope(rcycl.ppp, Lest, nsim = 99, rank = 1, global = T, correction = "Ripley")
plot(rLfun.global)

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

####DEPENDENCE####

##cdf null hypothesis - CDF of the covariate at all points is equal 
## to the CDF of covariate evaluated at the location of the point pattern
get_dependence_results = function(rcycl.ppp, covar, covar_string, test_string) {
  #plot(cdf.test(rcycl.ppp, covar))
  if(str_detect(test_string, "less")) {
    return(c(
      covar_string, test_string, 
      cdf.test(rcycl.ppp, covar, alternative = "less")['p.value'], 
      auc(rcycl.ppp, covar)
    ))
  } else if(str_detect(test_string, "greater")) {
    return(c(
      covar_string, test_string, 
      cdf.test(rcycl.ppp, covar, alternative = "greater")['p.value'], 
      auc(rcycl.ppp, covar)
    ))
  }else {
    return(c(
      covar_string, test_string, 
      cdf.test(rcycl.ppp, covar)['p.value'], 
      auc(rcycl.ppp, covar)
    ))
    
  }
}


dr = data.frame(
  covariate = character(), 
  test = character(), 
  p.val = integer(), 
  auc = integer()
)

##### Dependence of recycling points intensity on elevation, slope, erosion risk #####
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, dem_im, "DEM", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, dem_im, "DEM", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, dem_im, "DEM", "one-sided: greater")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sp_im, "slope percentage", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sp_im, "slope percentage", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sp_im, "slope percentage", "one-sided: greater")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sd_im, "slope degrees", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sd_im, "slope degrees", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sd_im, "slope degrees", "one-sided: greater")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, er_im, "erosion risk", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, er_im, "erosion risk", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, er_im, "erosion risk", "one-sided: greater")

##### Dependence of recycling points intensity on underlying artifact density ####
D = density(ppp, sigma=bw.diggle)
plot(D, useRaster=F)

dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(D), "artifact density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(D), "artifact density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(D), "artifact density", "one-sided: greater")


##### Dependence of recycling points intensity on underlying density of retouched artifacts ####
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(retouch.dens), "retouched artifact density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(retouch.dens), "retouched artifact density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(retouch.dens), "retouched artifact density", "one-sided: greater")

##### Volumetric and weight dependence #####
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(weight.dens), "artifact weight density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(weight.dens), "artifact weight density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(weight.dens), "artifact weight density", "one-sided: greater")

dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(thick.dens), "artifact thickness density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(thick.dens), "artifact thickness density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(thick.dens), "artifact thickness density", "one-sided: greater")

dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(length.dens), "artifact length density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(length.dens), "artifact length density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(length.dens), "artifact length density", "one-sided: greater")

dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(width.dens), "artifact width density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(width.dens), "artifact width density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(width.dens), "artifact width density", "one-sided: greater")

##### Artifact type dependence #####
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(compl_flk.dens), "complete flake density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(compl_flk.dens), "complete flake density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(compl_flk.dens), "complete flake density", "one-sided: greater")

dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(broke_flk.dens), "broken flake density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(broke_flk.dens), "broken flake density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(broke_flk.dens), "broken flake density", "one-sided: greater")

dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(core.dens), "core density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(core.dens), "core density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(core.dens), "core density", "one-sided: greater")

dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(tool.dens), "tool density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(tool.dens), "tool density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(tool.dens), "tool density", "one-sided: greater")


##### Weathering class dependence #####
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(str_weather.dens), "strongly weathered artifact density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(str_weather.dens), "strongly weathered artifact density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(str_weather.dens), "strongly weathered artifact density", "one-sided: greater")

dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(mid_weather.dens), "mildly weathered artifact density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(mid_weather.dens), "mildly weathered artifact density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(mid_weather.dens), "mildly weathered artifact density", "one-sided: greater")

dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(weak_weather.dens), "weakly weathered artifact density", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(weak_weather.dens), "weakly weathered artifact density", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, as.im(weak_weather.dens), "weakly weathered artifact density", "one-sided: greater")

dr$p.adj = p.adjust(dr$p.val, method = "BH")
dr$signif = ifelse(dr$p.adj < 0.05, TRUE, FALSE)

dr$covariate = factor(dr$covariate, 
                      levels = c(
                        "artifact density", "retouched artifact density",
                        "complete flake density", "broken flake density", 
                        "tool density", "core density", 
                        "strongly weathered artifact density", "mildly weathered artifact density", "weakly weathered artifact density", 
                        "artifact length density", "artifact width density", "artifact thickness density", "artifact weight density", 
                        "DEM", "slope percentage", "slope degrees", "erosion risk"
                      ))

dplot = ggplot(dr %>% filter(signif == TRUE) %>% filter(test != "two-sided")) +
  geom_col(aes(x = covariate, y = auc)) +
  geom_abline(aes(intercept = 0.5, slope = 0), color = "red") +
  coord_flip() +
  facet_wrap(~test, nrow = 1)
plot(dplot)
ggsave(filename = "figures/p2-dependence-auc.tiff", dpi = 300)

#### Marked point process for segregation analysis ####
marks(rcycl.ppp) = rcycl.data$Artfct_t
plot(rcycl.ppp)
plot(unmark(rcycl.ppp))
summary(rcycl.ppp)
plot(density(split(rcycl.ppp)), useRaster = F)

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


#### Figure -- recycling types spatial ####
sp.p2 = readOGR("data/artifact-shapefiles", "p2-artifacts")

sp.p2@data = sp.p2@data %>%
  mutate(Recycling_type =
           ifelse(str_detect(Rcyclng_n, "none"), "none", 
                  ifelse(str_detect(Rcyclng_n, " "), "multiple", 
                         Rcyclng_n))) %>%
  mutate(Recycling_type = str_replace_all(Recycling_type, "_", " ")) %>%
  mutate(Recycling_type = ifelse(Recycling_type == "core on flake blade", "core on flake/blade", Recycling_type))
sp.p2@data$Recycling_type = factor(sp.p2@data$Recycling_type, 
                                   levels = c("double patina", "core on flake/blade", "core on hammerstone", "multiple", "other", "none"))

sf.p2 = st_transform(st_as_sf(sp.p2), 32642)

pal1 = c("double patina" = "#0072B2", 
         "core on flake/blade" = "#D55E00", 
         "core on hammerstone" = "#CC79A7", 
         "multiple" = "#009E73", 
         "other" = "#F0E442")

plot.sf.p2 = sf.p2[!is.na(sf.p2$Recycling_type),]

plot.dem = as(dem_crop, "SpatialPixelsDataFrame")
dem.df = as.data.frame(plot.dem)
colnames(dem.df) = c("value", "x", "y")

plot.slope = as(sd_crop, "SpatialPixelsDataFrame")
slope.df = as.data.frame(plot.slope)
colnames(slope.df) = c("value", "x", "y")

rt.spat.plot = ggplot() +
  geom_tile(data = dem.df, aes(x = x, y = y, fill = value), alpha = 0.25) +
  geom_sf(data = plot.sf.p2[plot.sf.p2$Recycling_type == "none",], size = 0.25, alpha = 0.25) +
  geom_sf(data = plot.sf.p2[plot.sf.p2$Recycling_type != "none",], aes(color = Recycling_type), size = 1) +
  scale_fill_gradientn(colors = terrain.colors(10)) +
  scale_color_manual(values = pal1) +
  labs(color = "recycling signature", fill = "elevation", x = "", y = "") +
  theme(axis.text = element_blank()) +
  annotation_scale(pad_x = unit(0.75, "cm"), pad_y = unit(0.75, "cm"))
plot(rt.spat.plot)
ggsave(plot = rt.spat.plot, filename = "figures/p2/p2-recycling-signature-spatial_elev.tiff", dpi = 300)


rt.spat.plot2 = ggplot() +
  geom_tile(data = slope.df, aes(x = x, y = y, fill = value), alpha = 0.25) +
  geom_sf(data = plot.sf.p2[plot.sf.p2$Recycling_type == "none",], size = 0.25, alpha = 0.25) +
  geom_sf(data = plot.sf.p2[plot.sf.p2$Recycling_type != "none",], aes(color = Recycling_type), size = 1) +
  scale_fill_gradientn(colors = cm.colors(10)) +
  scale_color_manual(values = pal1) +
  labs(color = "recycling signature", fill = "slope (degrees)", x = "", y = "") +
  theme(axis.text = element_blank()) +
  annotation_scale(pad_x = unit(0.75, "cm"), pad_y = unit(0.75, "cm"))
plot(rt.spat.plot2)
ggsave(plot = rt.spat.plot2, filename = "figures/p2/p2-recycling-signature-spatial_slope.tiff", dpi = 300)

ggsave(filename = "figures/p2-map.tiff", ggarrange(rt.spat.plot, rt.spat.plot2), 
       dpi = 300, height = 15, width = 12)
