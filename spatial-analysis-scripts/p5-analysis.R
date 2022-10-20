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


Z = nnclean(rcycl.ppp, k=10, plothist = T)
plot(Z)

Kest = Kest(rcycl.ppp, correction = "best")
plot(Kest)

Lest = Lest(rcycl.ppp, correction = "best")
plot(Lest)
#both K function and L function indicate recycled objects are clustered

#### modeling poisson point processes ####
ppm0 = ppm(rcycl.ppp ~ 1)

ppm1 = ppm(rcycl.ppp ~ artifact.dens)
ppm2 = ppm(rcycl.ppp ~ log(artifact.dens))
ppm2

AIC(ppm1)
AIC(ppm2)
# power relationship fits better

anova(ppm0, ppm2, test = "LRT")
#artifact density significantly affects intensity of recycled artifacts

lambda0 = predict(ppm2)
rh1 = rhohat(rcycl.ppp, artifact.dens, baseline = lambda0)
plot(rh1)
#underestimates intensity of recycled artifacts at low artifact density areas,
#possibly slighly overestimates intensity of recycled artifacts at high artifact density areas

fit1 = ppm(rcycl.ppp ~ log(artifact.dens) +
             compl_flk.dens + broke_flk.dens +
             tool.dens + tool_frag.dens + 
             core.dens + core_frag.dens + shatter.dens)
fit1 #only tools and cores have significant effect, but with other things include effect of artifact density is not significant

anova(ppm2, fit1, test = "LRT") ## adding artifact types to the model DOES add significant information
drop1(fit1)
##adds information but none of the artifact types have a significant effect on the intensity of 
#recycled artifacts apart from artifact density effects

fit2 = ppm(rcycl.ppp ~ log(artifact.dens) +
      str_weather.dens + mid_weather.dens + weak_weather.dens + not_weather.dens)
fit2
anova(ppm2, fit2, test = "LRT") ## adding weathering to the model DOES add significant information
#significant effects of strongly weathered and weakly weathered artifacts

fit3 = ppm(rcycl.ppp ~ log(artifact.dens) +
      rside_dorsal.dens + rside_ventral.dens + rside_bifacial.dens)
fit3
anova(ppm2, fit3, test = "LRT") ## adding retouch side DOES add significant information
#significant effects of dorsal and bifacial retouch side counts

fit4 = ppm(rcycl.ppp ~ log(artifact.dens) + retouch.dens)
anova(ppm2, fit4, test = "LRT") ## adding retouch in general DOES add significant information
fit4 #significant effect of retouched artifacts

fit5 = ppm(rcycl.ppp ~ log(artifact.dens) + edge_dam.dens)
anova(ppm2, fit5, test = "LRT") ## adding edge damage DOES add significant information 
fit5 #significant effect

fit6 = ppm(rcycl.ppp ~ log(artifact.dens) +
             type_flake.dens + type_blade.dens + type_bladelet.dens)
anova(ppm2, fit6, test = "LRT") # blank type DOES add significant information
fit6 #significant effect of flake blank intensity

fit7 = ppm(rcycl.ppp ~ log(artifact.dens) +
             ttype_mult.dens + ttype_notch.dens + ttype_dent.dens +
             ttype_nd.dens + ttype_point.dens + ttype_oth.dens)
fit7
anova(ppm2, fit7, test = "LRT") #adding tool type does add significant information, but none of the tool types have individual significant effects

fit8 = ppm(rcycl.ppp ~ log(artifact.dens) + length.dens + width.dens + thick.dens + weight.dens)
anova(ppm2, fit8, test = "LRT") #adding all size covariates together does not add significant information
fit8
fit8.1 = ppm(rcycl.ppp ~ log(artifact.dens) + length.dens)
anova(ppm2, fit8.1, test = "LRT") #adding length individually adds significant information
fit8.1
fit8.2 = ppm(rcycl.ppp ~ log(artifact.dens) + width.dens)
anova(ppm2, fit8.2, test = "LRT") #adding width individually adds significant information
fit8.2
fit8.3 = ppm(rcycl.ppp ~ log(artifact.dens) + thick.dens)
anova(ppm2, fit8.3, test = "LRT") #adding thickness individually adds significant information
fit8.3
fit8.4 = ppm(rcycl.ppp ~ log(artifact.dens) + weight.dens)
anova(ppm2, fit8.4, test = "LRT") #adding weight individually adds significant information
fit8.4

fit9 = ppm(rcycl.ppp ~ log(artifact.dens) + cortex.dens)
anova(ppm2, fit9, test = "LRT") #adding cortex does not add significant information


full.ppm = ppm(rcycl.ppp ~ log(artifact.dens) + 
                 compl_flk.dens + broke_flk.dens +
                 tool.dens + tool_frag.dens + 
                 core.dens + core_frag.dens + shatter.dens +
                 str_weather.dens + mid_weather.dens + weak_weather.dens + not_weather.dens +
                 retouch.dens +
                 rside_dorsal.dens + rside_ventral.dens + rside_bifacial.dens +
                 edge_dam.dens +
                 type_flake.dens + type_blade.dens + type_bladelet.dens +
                 ttype_mult.dens + ttype_notch.dens + ttype_dent.dens +
                 ttype_nd.dens + ttype_point.dens + ttype_oth.dens )
full.ppm
