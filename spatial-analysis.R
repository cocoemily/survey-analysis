# spatial analysis
#point pattern data


library(sp)
library(rgdal)
library(spdep)
library(tidyverse)
library(sf)

#analysis of point pattern data
library(spatstat)
library(splancs)
library(maptools)

sp.artifacts = readOGR("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/", "june-artifacts")
sp.artifacts@data$recycled = !is.na(sp.artifacts@data$Rcyclng_d)
sp.artifacts@data = sp.artifacts@data %>% mutate(location = ifelse(St_nm %in% c("Square 1", "Square  1", "Square 2", "Square 3"), "P1", 
                                                  ifelse(St_nm %in% c("Square 4", "Square 5"), "P2", "P5")))

#plot(sp.artifacts)

sf.artifacts = st_as_sf(sp.artifacts)

ggplot() +
  geom_sf(data = sf.artifacts, aes(fill = recycled))


sp.p1 = readOGR("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/artifact-shapefiles", "P1-june")
sp.p2 = readOGR("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/artifact-shapefiles", "P2-june")
sp.p5 = readOGR("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/artifact-shapefiles", "P5-june")

sf.p1 = st_as_sf(sp.p1)
sf.p2 = st_as_sf(sp.p2)
sf.p5 = st_as_sf(sp.p5)

ggplot() +
  geom_sf(data = sf.p1, aes(color = as.factor(rcycl)))

ggplot() +
  geom_sf(data = sf.p2, aes(color = as.factor(rcycl)))

ggplot() +
  geom_sf(data = sf.p5, aes(color = as.factor(rcycl)))


xy = coordinates(sp.p1)
w = owin(xrange = c(min(xy[,1]), max(xy[,1])), yrange = c(min(xy[,2]), max(xy[,2])))


as.ppp(sf.p1)





kn = knearneigh(sf.artifacts, k = 3)
art.nb = knn2nb(kn)
art.lw = nb2listw(art.nb)

moran.test(sf.artifacts$Weght, art.lw)

#morans I on true/false data?
#morans I on categorical data?


I = localmoran(sf.artifacts$Weght, art.lw)



##cluster analysis




