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
source("spatial-analysis-scripts/get-covariates.R")

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))

rcycl.data = st_data %>% filter(rcycl == 1)
rcycl.data = rcycl.data %>% filter(poss_roll == F | is.na(poss_roll))
rcycl.ppp = as.ppp(st_data %>% filter(rcycl == 1))
marks(rcycl.ppp) = NULL
Window(rcycl.ppp) = win
plot(rcycl.ppp)

K2 = density(rcycl.ppp, sigma = bw.diggle, adjust = 2)
plot(K2, useRaster = F, main = "Recycled artifact density")
contour(K2, add = TRUE)


#### recycling point process and artifact density ####
cdf.test(rcycl.ppp, artifact.dens, test = "ks")
#p value less than 0.05, so reject the null -> recycled points do depend on underlying density of artifacts
berman.test(rcycl.ppp, artifact.dens)
#recycled points depend on underlying density of artifacts

coproc = roc(rcycl.ppp, artifact.dens)
plot(coproc)
#artifact density has very strong discriminatory power, strong effect of artifact density on recycled object point process

##nearest neighbor cleaning -- separates noise from features 
Z = nnclean(rcycl.ppp, k=10, plothist = T)
plot(Z)

Kest = Kest(rcycl.ppp, correction = "best")
plot(Kest)

Lest = Lest(rcycl.ppp, correction = "best")
plot(Lest)
#both K function and L function indicate recycled objects are more dispersed than would be expected

E = envelope(rcycl.ppp, fun = "Kest", nsim = npoints(rcycl.ppp), fix.n = T)
plot(E)


#### modeling poisson point processes ####
ppm0 = ppm(rcycl.ppp ~ 1)

ppm1 = ppm(rcycl.ppp ~ artifact.dens)
ppm2 = ppm(rcycl.ppp ~ log(artifact.dens))
ppm2

AIC(ppm1)
AIC(ppm2)

#power relationships between recycled object intensity and artifact density 

anova(ppm0, ppm2, test = "LRT")
#artifact density significantly affects intensity of recycled artifacts

lambda0 = predict(ppm2)
rh1 = rhohat(rcycl.ppp, artifact.dens, baseline = lambda0)
plot(rh1)
#underestimates intensity of recycled artifacts at low artifact density areas

fit1 = ppm(rcycl.ppp ~ log(artifact.dens) +
             compl_flk.dens + broke_flk.dens +
             tool.dens + tool_frag.dens + 
             core.dens)
fit1 #only tools and cores have significant effect, but with other things include effect of artifact density is not significant

anova(ppm2, fit1, test = "LRT") ## adding artifact types to the model DOES add significant information
drop1(fit1)

fit1.1 = ppm(rcycl.ppp ~ log(artifact.dens) + tool.dens + core.dens) 
anova(ppm2, fit1.1, test = "LRT")
fit1.1


fit3 = ppm(rcycl.ppp ~ log(artifact.dens) +
             str_weather.dens + mid_weather.dens + weak_weather.dens + not_weather.dens)
fit3 # only weakly weathered has a signficant effect
anova(ppm2, fit3, test = "LRT") ## adding weathering does not add significant information

fit4 = ppm(rcycl.ppp ~ log(artifact.dens) +
             rside_dorsal.dens + rside_ventral.dens + rside_bifacial.dens)
fit4
anova(ppm2, fit4, test = "LRT") ## adding retouch side does not add signficant information

fit5 = ppm(rcycl.ppp ~ log(artifact.dens) + retouch.dens)
fit5
anova(ppm2, fit5, test = "LRT") ## adding retouched artifact density does not add significant information

fit6 = ppm(rcycl.ppp ~ log(artifact.dens) + edge_dam.dens)
fit6
anova(ppm2, fit6, test = "LRT") ## adding edge damaged artifact density does not add significant information

fit7 = ppm(rcycl.ppp ~ log(artifact.dens) +
             type_flake.dens + type_blade.dens + type_bladelet.dens)
fit7
anova(ppm2, fit7, test = "LRT") ## adding blank type does not add significant information

fit8 = ppm(rcycl.ppp ~ log(artifact.dens) +
             ttype_mult.dens + ttype_notch.dens + ttype_dent.dens + ttype_oth.dens)
fit8
anova(ppm2, fit8, test = "LRT") ## adding tool type does not add significant information

fit9 = ppm(rcycl.ppp ~ log(artifact.dens) + length.dens + width.dens + thick.dens + weight.dens)
fit9
drop1(fit9)
fit9.1 = ppm(rcycl.ppp ~ log(artifact.dens) + weigth.dens)
fit9.1
anova(ppm2, fit9.1, test = "LRT")
anova(ppm2, fit9, test = "LRT")
## adding in size information does not add significant information


