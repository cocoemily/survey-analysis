## Poisson point processes

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

set.seed(120109)

p1.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p1-artifacts")
p2.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p2-artifacts")
p5.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p5-artifacts")

p1.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p1-sqs")
p2.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p2-sqs")
p5.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p5-sqs")


data = p1.artifacts
window = p1.window

source("spatial-analysis-scripts/clean-data.R")
source("spatial-analysis-scripts/get-covariates.R")


#### SPATIAL POINT PATTERN ANALYSIS IN R TESTING ####
st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))


#### POINT PROCESS ANALYSIS ####
#### ppp with recycled artifacts ####
rcycl.ppp = as.ppp(st_data %>% filter(rcycl == 1))
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
LR = scanLRTS(rcycl.ppp, r = 2 * bw.diggle(rcycl.ppp))
plot(LR, useRaster = F)
pvals = eval.im(pchisq(LR, df = 1, lower.tail = F))
plot(pvals, useRaster = F)
plot(rcycl.ppp, add = T, col = "white")

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
#results show that estimated intensity of recycled artifacts when artifact density is 0:
#0.045 per square unit
#intensity of recycled artifacts increases by a factor of 4 if artifact density increase to 1
plot(effectfun(ppm1, "artifact.dens", se.fit = T))
#with this model the dependence is tightly specified 
#--> assumes intensity of recycled artifacts is expotentially related to artifact density

#power law relationships between recycled artifacts and artifact density?
ppm2 = ppm(rcycl.ppp ~ log(artifact.dens))
ppm2

AIC(ppm1)
AIC(ppm2)
#AIC criterion suggests that ppm1 is better with power law relationship between values


anova(ppm0, ppm2, test = "LR")
#artifact density significantly affects intensity of recycled artifacts

diagnose.ppm(ppm0)


lam0 = fitted(ppm2, dataonly = T)
rcycl.dens = density(rcycl.ppp, weights = 1/lam0)
range(rcycl.dens)
#wide range of relative intensity --> poor fit of ppm1

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


fit3 = ppm(rcycl.ppp ~ log(artifact.dens) +
             str_weather.dens + mid_weather.dens + weak_weather.dens + not_weather.dens)
fit3
AIC(fit3)
drop1(fit3) #model performs best without mildly weathered covariate

fit4 = ppm(rcycl.ppp ~ log(artifact.dens) +
             rside_dorsal.dens + rside_ventral.dens + rside_bifacial.dens)
fit4
AIC(fit4)
drop1(fit4) #including retouch side does perform better than model just based on artifact density

fit5 = ppm(rcycl.ppp ~ log(artifact.dens) +
             type_flake.dens + type_blade.dens + type_bladelet.dens)
fit5
AIC(fit5)
drop1(fit5) #include blank type does not perform better than base model

fit6 = ppm(rcycl.ppp ~ log(artifact.dens) +
             ttype_mult.dens + ttype_notch.dens + ttype_dent.dens +
             ttype_nd.dens + ttype_point.dens + ttype_oth.dens)
fit6
AIC(fit6)
drop1(fit6) #only need to keep multiple tool types


p1.ffit = ppm(rcycl.ppp ~ log(artifact.dens) +
             compl_flk.dens + broke_flk.dens +
             tool.dens + tool_frag.dens + 
             str_weather.dens +  weak_weather.dens + not_weather.dens +
             ttype_mult.dens)
p1.ffit
AIC(p1.ffit)
