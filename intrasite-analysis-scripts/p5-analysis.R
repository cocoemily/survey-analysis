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

theme_set(theme_bw())

p5.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p5-artifacts")
p5.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p5-sqs")

data = p5.artifacts
window = p5.window

source("intrasite-analysis-scripts/clean-data.R")
#source("intrasite-analysis-scripts/get-covariates.R")

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))

####ALL POINTS####
quadrat.test(ppp, 5, method = "MonteCarlo") 
#results of quadrat test indicate a inhomogeneous process
plot(envelope(ppp, fun = Ginhom, nsim = 99))
plot(envelope(ppp, fun = Finhom, nsim = 99))
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
# Lin = Linhom(ppp, lambda = Dnn)
# plot(Lin)
#both K function falls slightly below the poisson line --> points are slightly more dispersed than expected


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
wdata = st_data %>% filter(!is.na(Wthrng_c))
ppp = as.ppp(wdata)
marks(ppp) = wdata$Wthrng_c
Window(ppp) = win
plot(ppp)
segregation.test(ppp, nsim=99) 
#there is spatial segregation of weathering types
plot(density(split(ppp)), useRaster=F)

# #weathering class without other weathering
# wdata = st_data %>% filter(Wthrng_c != "other" & !is.na(Wthrng_c))
# wppp = as.ppp(wdata)
# marks(wppp) = wdata$Wthrng_c
# Window(wppp) = win
# plot(wppp)
# segregation.test(wppp, nsim=99)
# plot(density(split(wppp)), useRaster=F)
# #spatial segregation of weathering types

ttdata = st_data %>% filter(!is.na(tool.type))
ttppp = as.ppp(ttdata)
marks(ttppp) = ttdata$tool.type
Window(ttppp) = win
plot(ttppp)
segregation.test(ttppp, nsim=99)
#spatial segregation of tool types
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


####DEPENDENCE####

##cdf null hypothesis - CDF of the covariate at all points is equal 
## to the CDF of covariate evaluated at the location of the point pattern
get_dependence_results = function(rcycl.ppp, covar, covar_string, test_string) {
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

##### Dependence of recycling points intensity on elevation and slope #####
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, dem_im, "DEM", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, dem_im, "DEM", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, dem_im, "DEM", "one-sided: greater")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sp_im, "slope percentage", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sp_im, "slope percentage", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sp_im, "slope percentage", "one-sided: greater")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sd_im, "slope degrees", "two-sided")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sd_im, "slope degrees", "one-sided: less")
dr[nrow(dr) + 1, ] <- get_dependence_results(rcycl.ppp, sd_im, "slope degrees", "one-sided: greater")

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

dr$signif = ifelse(dr$p.val < 0.05, TRUE, FALSE)
dr$covariate = factor(dr$covariate, 
                      levels = c(
                        "artifact density", "retouched artifact density",
                        "complete flake density", "broken flake density", 
                        "tool density", "core density", 
                        "strongly weathered artifact density", "mildly weathered artifact density", "weakly weathered artifact density", 
                        "artifact length density", "artifact width density", "artifact thickness density", "artifact weight density", 
                        "DEM", "slope percentage", "slope degrees"
                      ))

dplot = ggplot(dr %>% filter(signif == TRUE) %>% filter(test != "two-sided")) +
  geom_col(aes(x = covariate, y = auc)) +
  geom_abline(aes(intercept = 0.5, slope = 0), color = "red") +
  coord_flip() +
  facet_wrap(~test, nrow = 1)
plot(dplot)
ggsave(filename = "figures/p5-dependence-auc.tiff", dpi = 300)


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

sp.p5@data = sp.p5@data %>%
  mutate(Recycling_type =
           ifelse(str_detect(Rcyclng_n, "none"), "none", 
                  ifelse(str_detect(Rcyclng_n, " "), "multiple", 
                         Rcyclng_n))) %>%
  mutate(Recycling_type = str_replace_all(Recycling_type, "_", " ")) %>%
  mutate(Recycling_type = ifelse(Recycling_type == "core on flake blade", "core on flake/blade", Recycling_type))
sp.p5@data$Recycling_type = factor(sp.p5@data$Recycling_type, 
                                   levels = c("double patina", "core on flake/blade", "core on hammerstone", "multiple", "other", "none"))

sf.p5 = st_transform(st_as_sf(sp.p5), 32642)

pal1 = c("double patina" = "#0072B2", 
         "core on flake/blade" = "#D55E00", 
         "core on hammerstone" = "#CC79A7", 
         "multiple" = "#009E73", 
         "other" = "#F0E442")

plot.sf.p5 = sf.p5[!is.na(sf.p5$Recycling_type),]

plot.dem = as(dem_crop, "SpatialPixelsDataFrame")
dem.df = as.data.frame(plot.dem)
colnames(dem.df) = c("value", "x", "y")

plot.slope = as(sd_crop, "SpatialPixelsDataFrame")
slope.df = as.data.frame(plot.slope)
colnames(slope.df) = c("value", "x", "y")

rt.spat.plot = ggplot() +
  geom_tile(data = dem.df, aes(x = x, y = y, fill = value), alpha = 0.25) +
  geom_sf(data = plot.sf.p5[plot.sf.p5$Recycling_type == "none",], size = 0.25, alpha = 0.25) +
  geom_sf(data = plot.sf.p5[plot.sf.p5$Recycling_type != "none",], aes(color = Recycling_type), size = 1) +
  scale_fill_gradientn(colors = terrain.colors(10)) +
  scale_color_manual(values = pal1) +
  labs(color = "recycling signature", fill = "elevation", x = "", y = "") +
  theme(axis.text = element_blank())
ggsave(plot = rt.spat.plot, filename = "figures/p5-recycling-signature-spatial_elev.tiff", dpi = 300)

plot(rt.spat.plot)

rt.spat.plot2 = ggplot() +
  geom_tile(data = slope.df, aes(x = x, y = y, fill = value), alpha = 0.25) +
  geom_sf(data = plot.sf.p5[plot.sf.p5$Recycling_type == "none",], size = 0.25, alpha = 0.25) +
  geom_sf(data = plot.sf.p5[plot.sf.p5$Recycling_type != "none",], aes(color = Recycling_type), size = 1) +
  scale_fill_gradientn(colors = cm.colors(10)) +
  scale_color_manual(values = pal1) +
  labs(color = "recycling signature", fill = "slope (degrees)", x = "", y = "") +
  theme(axis.text = element_blank())
plot(rt.spat.plot2)
ggsave(plot = rt.spat.plot2, filename = "figures/p5-recycling-signature-spatial_slope.tiff", dpi = 300)



