# spatial analysis
#point pattern data

library(sp)
library(rgdal)
library(spdep)
library(tidyverse)
library(sf)
library(raster)

#analysis of point pattern data
library(spatstat)
library(splancs)
library(maptools)

theme_set(theme_bw())

# sp.artifacts = readOGR("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/", "june-artifacts")
# sp.artifacts@data$recycled = !is.na(sp.artifacts@data$Rcyclng_d)
# sp.artifacts@data = sp.artifacts@data %>% mutate(location = ifelse(St_nm %in% c("Square 1", "Square  1", "Square 2", "Square 3"), "P1", 
#                                                   ifelse(St_nm %in% c("Square 4", "Square 5"), "P2", "P5")))
# 
# #plot(sp.artifacts)
# 
# sf.artifacts = st_as_sf(sp.artifacts)
# 
# ggplot() +
#   geom_sf(data = sf.artifacts, aes(color = recycled))


sp.p1 = readOGR("data/artifact-shapefiles", "p1-artifacts2")
p1.counts = readOGR("data/practice", "recycling_counts")
#sp.p2 = readOGR("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/artifact-shapefiles", "P2-june")
#sp.p5 = readOGR("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/artifact-shapefiles", "P5-june")


sf.p1 = st_as_sf(sp.p1)
#sf.p2 = st_as_sf(sp.p2)
#sf.p5 = st_as_sf(sp.p5)

ggplot() +
  geom_sf(data = sf.p1, aes(color = Rcyclng_n))

# ggplot() +
#   geom_sf(data = sf.p2, aes(color = as.factor(rcycl)))
# 
# ggplot() +
#   geom_sf(data = sf.p5, aes(color = as.factor(rcycl)))


xy = coordinates(sp.p1)
win.sp = readOGR("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/artifact-shapefiles", "p1-window")
win = st_read("~/Desktop/NYU/Dissertation-Research/Survey/June-survey/artifact-shapefiles", "p1-window")
win2 = st_transform(win, crs = 6345)
w = as.owin(win2)
p1_proj = st_transform(sf.p1, crs = 6345)
xy = st_coordinates(p1_proj)

ppp = ppp(x = xy[,1], y = xy[,2], window = w)

plot(ppp)

ds = density(ppp)
plot(ds)


##issue with the extent shapefile for square 6 - some artifacts are outside of the extent



kn = knearneigh(sp.p1, k = 3)
art.nb = knn2nb(kn)
art.lw = nb2listw(art.nb)

moran.test(sp.p1$Weght, art.lw)
globalG.test(sp.p1$Weght,art.lw,alternative="two.sided")

#morans I on true/false data?
#morans I on categorical data?


I = localmoran(sp.p1$Weght, art.lw)
G = localG(sp.p1$Weght, art.lw)



local_g.ma=as.matrix(G)
sp.p1 = cbind(sf.p1,local_g.ma)

sp.p1$pval =  2*pnorm(-abs(sp.p1$local_g.ma))
subset(sp.p1,local_g.ma>0&pval<0.05)$Id_nm
sp.p1$sign = ifelse(sp.p1$pval<0.05, "p<0.05", "p>0.05")

ggplot(sp.p1) +
  geom_sf(aes(color = local_g.ma, shape = sign)) +
  scale_color_viridis_c()

#hotspot analysis

# pixelsize = 100
# box = round(extent(sf.p1) / pixelsize) * pixelsize
# template = raster(box, crs = 6345,
#                   nrows = (box@ymax - box@ymin) / pixelsize, 
#                   ncols = (box@xmax - box@xmin) / pixelsize)
  

ext = extent(p1.counts)
plot(p1.counts)
temp = raster(p1.counts, ext = ext)
test = rasterize(p1.counts, temp, field = "NUMPOINTS")

plot(test)


##cluster analysis



#exploratory plots
sp.p1 = sp.p1 %>%
  mutate(Recycling_type =
           ifelse(str_detect(Rcyclng_n, "none"), "none", 
                  ifelse(str_detect(Rcyclng_n, " "), "multiple", 
                         Rcyclng_n))) %>%
  mutate(Recycling_type = str_replace_all(Recycling_type, "_", " ")) %>%
  mutate(Recycling_type = ifelse(Recycling_type == "core on flake blade", "core on flake/blade", Recycling_type))
sp.p1$Recycling_type = factor(sp.p1$Recycling_type, 
                                   levels = c("double patina", "core on flake/blade", "core on hammerstone", "multiple", "other", "none"))

lower_bound = median(sp.p1$Flk_l, na.rm = T) - 3 * mad(sp.p1$Flk_l, constant = 1, na.rm = T)
upper_bound = median(sp.p1$Flk_l, na.rm = T) + 3 * mad(sp.p1$Flk_l, constant = 1, na.rm = T)


flk.p1 = sp.p1 %>% filter(Flk_l <= upper_bound & Flk_l >= lower_bound)
ggplot(flk.p1) +
  geom_sf(aes(color = Flk_l, shape = Recycling_type)) +
  scale_color_viridis_c()

ggplot(data = p1.counts, aes(fill = NUMPOINTS)) +
  geom_polygon()



