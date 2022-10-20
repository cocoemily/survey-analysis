## Poisson point process analysis for Semizbugu P1
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

p1.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p1-artifacts")
p1.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p1-sqs")

data = p1.artifacts
window = p1.window

source("spatial-analysis-scripts/clean-data.R")
source("spatial-analysis-scripts/get-covariates.R")

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))


#### POINT PROCESS ANALYSIS ####
#### ppp with recycled artifacts ####
rcycl.data = st_data %>% filter(rcycl == 1)
rcycl.data = rcycl.data %>% filter(poss_roll == F | is.na(poss_roll))
rcycl.ppp = as.ppp(rcycl.data)
marks(rcycl.ppp) = NULL
Window(rcycl.ppp) = win
plot(rcycl.ppp)

K2 = density(rcycl.ppp, sigma = bw.diggle, adjust = 2)
plot(K2, useRaster = F, main = "Recycled artifact density")
contour(K2, add = TRUE)

#test for complete spatial randomness
plot(Gest(rcycl.ppp))

plot(Fest(rcycl.ppp))
plot(Jest(rcycl.ppp))
#all measures confirm CSR --> homogenous process

#### recycling point process and artifact density ####
cdf.test(rcycl.ppp, artifact.dens, test = "ks")
#p value less than 0.05, so reject the null -> recycled points do depend on underlying density of artifacts
berman.test(rcycl.ppp, artifact.dens)
#recycled points depend on underlying density of artifacts


coproc = roc(rcycl.ppp, artifact.dens)
plot(coproc)
#artifact density has very strong discriminatory power, strong effect of artifact density on recycled object point process


#find hotspots
# LR = scanLRTS(rcycl.ppp, r = 2 * bw.diggle(rcycl.ppp))
# plot(LR, useRaster = F)
# pvals = eval.im(pchisq(LR, df = 1, lower.tail = F))
# plot(pvals, useRaster = F)
# plot(rcycl.ppp, add = T, col = "white")

##nearest neighbor cleaning -- separates noise from features 
Z = nnclean(rcycl.ppp, k=10, plothist = T)
plot(Z)


Kest = Kest(rcycl.ppp, correction = "best")
plot(Kest)

Lest = Lest(rcycl.ppp, correction = "best")
plot(Lest)
#both K function and L function indicate a regular point process for recycled objects

E = envelope(rcycl.ppp, fun = "Kest", nsim = npoints(rcycl.ppp), fix.n = T)
plot(E)



#### modeling poisson point processes ####
ppm0 = ppm(rcycl.ppp ~ 1)

ppm1 = ppm(rcycl.ppp ~ artifact.dens)
ppm1

#power law relationships between recycled artifacts and artifact density?
ppm2 = ppm(rcycl.ppp ~ log(artifact.dens))
ppm2
plot(effectfun(ppm2, "artifact.dens", se.fit = T))

AIC(ppm1)
AIC(ppm2)
#AIC criterion suggests that ppm1 is better with power law relationship between values


anova(ppm0, ppm2, test = "LRT")
#artifact density significantly affects intensity of recycled artifacts


# lam0 = fitted(ppm2, dataonly = T)
# rcycl.dens = density(rcycl.ppp, weights = 1/lam0)
# range(rcycl.dens)
# #wide range of relative intensity --> poor fit of ppm2

lambda0 = predict(ppm2)
rh1 = rhohat(rcycl.ppp, artifact.dens, baseline = lambda0)
plot(rh1)
#underestimates intensity of recycled artifacts at highest artifact density
#but otherwise performs pretty well

fit1 = ppm(rcycl.ppp ~ log(artifact.dens) +
             log(compl_flk.dens) + log(broke_flk.dens) +
             log(tool.dens) + log(tool_frag.dens) + 
             log(core.dens) + log(core_frag.dens))
fit1

fit2 = ppm(rcycl.ppp ~ log(artifact.dens) +
             compl_flk.dens + broke_flk.dens +
             tool.dens + tool_frag.dens + 
             core.dens + core_frag.dens)
fit2

AIC(fit1)
AIC(fit2)
#do not log the values of separate artifact types

drop1(fit2) #model performs better without core and core fragment covariates

fit2.2 = ppm(rcycl.ppp ~ log(artifact.dens) + 
               compl_flk.dens + broke_flk.dens +
               tool.dens + tool_frag.dens )
fit2.2
drop1(fit2.2) #new model performs better with complete flake covariate

fit2.3 = ppm(rcycl.ppp ~ log(artifact.dens) + 
               broke_flk.dens +
               tool.dens + tool_frag.dens )
fit2.3
drop1(fit2.3)
AIC(fit2.3)

anova(ppm2, fit2, test = "LRT") ## adding in artifact types adds significant information
fit2


fit3 = ppm(rcycl.ppp ~ log(artifact.dens) +
             str_weather.dens + mid_weather.dens + weak_weather.dens + not_weather.dens)
fit3
AIC(fit3)
drop1(fit3) #model performs best without mildly weathered covariate

anova(ppm2, fit3, test = "LRT") ## adding weathering does not add significant information


fit4 = ppm(rcycl.ppp ~ log(artifact.dens) +
             rside_dorsal.dens + rside_ventral.dens + rside_bifacial.dens)
fit4
AIC(fit4)
drop1(fit4) #including retouch side does perform better than model just based on artifact density

anova(ppm2, fit4, test = "LRT")## adding retouch side does not add significant information

fit4.1 = ppm(rcycl.ppp ~ log(artifact.dens) + retouch.dens)
fit4.1
anova(ppm2, fit4.1, test = "LRT") ## adding retouch intensity does not add significant information


fit5 = ppm(rcycl.ppp ~ log(artifact.dens) +
             type_flake.dens + type_blade.dens + type_bladelet.dens)
fit5
AIC(fit5)
drop1(fit5) #include blank type does not perform better than base model

anova(ppm2, fit5, test = "LRT") ## adding blank type does not add significant information

fit6 = ppm(rcycl.ppp ~ log(artifact.dens) +
             ttype_mult.dens + ttype_notch.dens + ttype_dent.dens +
             ttype_nd.dens + ttype_point.dens + ttype_oth.dens)
fit6
AIC(fit6)
drop1(fit6) #only need to keep multiple tool types

anova(ppm2, fit6, test = "LRT") ## adding tool type does not add significant information


fit7 = ppm(rcycl.ppp ~ log(artifact.dens) + length.dens + width.dens + thick.dens + weight.dens)
fit7
AIC(fit7)
drop1(fit7) #none of the size covariates impact intensity of recycled artifacts 
anova(ppm2, fit7, test = "LRT") ## adding size covariates does not add significant information


fit8 = ppm(rcycl.ppp ~ log(artifact.dens) + cortex.dens)
anova(ppm2, fit8, test = "LRT") # cortex covariate does not add significant information

fit9 = ppm(rcycl.ppp ~ log(artifact.dens) + edge_dam.dens)
anova(ppm2, fit9, test = "LRT") # edge damage covariate does not add significant information


