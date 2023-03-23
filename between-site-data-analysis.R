library(tidyverse)
library(sp)
library(raster)
library(rgeos)
library(pscl)
library(spatstat)
library(maptools)
library(rcompanion)
library(rgdal)
library(sf)
library(ggspatial)

theme_set(theme_bw())

data = read_csv("data/between-site-points.csv")
reach_view_points <- read_csv("data/reach-view-points.csv")

aggs = data %>% left_join(reach_view_points, by = c("reach-view-num" = "Name"))

aggs$Date = strptime(aggs$Date, format = "%Y-%m-%d %H:%M:%S")
aggs$start = strptime(aggs$start, format = "%Y-%m-%d %H:%M:%S")

aggs$"_point_altitude" = aggs$`Ellipsoidal height`
aggs$"_point_longitude" = aggs$Longitude
aggs$"_point_latitude" = aggs$Latitude

aggs$point = paste(aggs$Latitude, aggs$Longitude, aggs$`Ellipsoidal height`)

write_csv(aggs, "/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/July-survey/cleaned_between-site-points.csv")

p1.area = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p1-sqs")
p1.rpj = spTransform(p1.area, CRS("+init=epsg:32642"))
p2.area = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p2-sqs")
p2.rpj = spTransform(p2.area, CRS("+init=epsg:32642"))
p5.area = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p5-sqs")
p5.rpj = spTransform(p5.area, CRS("+init=epsg:32642"))

xy = aggs[,c(7,6)]
spaggs = SpatialPointsDataFrame(coords = xy, 
                                data = aggs, 
                                proj4string = crs(p1.area))
plot(spaggs)
spaggs.rpj = spTransform(spaggs, CRS("+init=epsg:32642"))

spaggs.rpj$dist_to_p1 = t(gDistance(spaggs.rpj, p1.rpj, byid = T))[,1]
spaggs.rpj$dist_to_p2 = t(gDistance(spaggs.rpj, p2.rpj, byid = T))[,1]
spaggs.rpj$dist_to_p5 = t(gDistance(spaggs.rpj, p5.rpj, byid = T))[,1]

spaggs.rpj$recycled.objects = ifelse(is.na(spaggs.rpj$Number_recycled_objects), 0, 
                                     spaggs.rpj$Number_recycled_objects)
hist(spaggs.rpj$recycled.objects)

plot.aggs = st_transform(st_as_sf(spaggs.rpj), 32642)

b_window = st_buffer(plot.aggs, 100)

dem = raster("/Users/emilycoco/Desktop/NYU/Dissertation-Research/sat-imagery/Kazakhstan/terrain/WorldDEMNeo_DSM_015_N47_09_E077_75/DEM/WorldDEMNeo_DSM_015_N47_09_E077_75_DEM.tif")
dem_rpj = projectRaster(dem, crs = 32642)
dem_crop = crop(dem_rpj, b_window)
plot(dem_crop)
plot.dem = as(dem_crop, "SpatialPixelsDataFrame")
dem.df = as.data.frame(plot.dem)
colnames(dem.df) = c("value", "x", "y")

slope.deg = raster("/Users/emilycoco/Desktop/NYU/Dissertation-Research/sat-imagery/Kazakhstan/derived/slope-percentage.tif")
sd_rpj = projectRaster(slope.deg, crs = 32642)
sd_crop = crop(sd_rpj, b_window)
plot.slope = as(sd_crop, "SpatialPixelsDataFrame")
slope.df = as.data.frame(plot.slope)
colnames(slope.df) = c("value", "x", "y")

plot1 = ggplot() +
  geom_tile(data = dem.df, aes(x = x, y = y, fill = value), alpha = 0.5) +
  geom_sf(data =  st_transform(st_as_sf(p1.area), 32642)) +
  geom_sf(data = st_transform(st_as_sf(p2.area), 32642)) +
  geom_sf(data = st_transform(st_as_sf(p5.area), 32642)) +
  geom_sf(data = plot.aggs, aes(color = recycled.objects)) + 
  scale_colour_gradient(low = "white", high = "black") +
  scale_fill_gradientn(colors = terrain.colors(50)) +
  labs(color = "recycled artifact count", fill = "elevation", x = "", y = "") +
  theme(axis.text = element_text(size = 5), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.4, "cm")) +
  annotation_scale(pad_x = unit(1, "cm"), pad_y = unit(0.5, "cm"))

plot(plot1)
ggsave(filename = "figures/between-site-plot_elev.tiff", plot1, 
       dpi = 300, width = 8, height = 6)

plot2 = ggplot() +
  geom_tile(data = slope.df, aes(x = x, y = y, fill = value), alpha = 0.5) +
  geom_sf(data =  st_transform(st_as_sf(p1.area), 32642)) +
  geom_sf(data = st_transform(st_as_sf(p2.area), 32642)) +
  geom_sf(data = st_transform(st_as_sf(p5.area), 32642)) +
  geom_sf(data = plot.aggs, aes(color = recycled.objects)) + 
  scale_colour_gradient(low = "thistle", high = "black") +
  scale_fill_gradientn(colors = cm.colors(50)) +
  labs(color = "recycled artifact count", fill = "slope (degrees)", x = "", y = "") +
  theme(axis.text = element_text(size = 5), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.4, "cm")) +
  annotation_scale(pad_x = unit(1, "cm"), pad_y = unit(0.5, "cm"))
plot(plot2)

ggsave(filename = "figures/between-site-plot_slope.tiff", plot2, 
       dpi = 300, width = 8, height = 6)


spaggs.rpj$dem = extract(dem_rpj, spaggs.rpj)
spaggs.rpj$slope_deg = extract(sd_rpj, spaggs.rpj)

dist.data = spaggs.rpj@data %>%
  dplyr::select(recycled.objects, 
                dem, 
                slope_deg, 
                dist_to_p1, 
                dist_to_p2, 
                dist_to_p5)
dist.data$dist_to_sa = apply(dist.data[,c(4:6)], 1, FUN = min)


fit1 = glm(recycled.objects ~ dist_to_p1 + dist_to_p2 + dist_to_p5, data = dist.data, family = "poisson")
summary(fit1)
pR2(fit1)

fit2 = glm(recycled.objects ~ dem, data = dist.data, family = "poisson")
summary(fit2)
pR2(fit2)

fit3 = glm(recycled.objects ~ slope_deg, data = dist.data, family = "poisson")
summary(fit3)
pR2(fit3)

fit4 = glm(recycled.objects ~ dist_to_sa, data = dist.data, family = "poisson")
summary(fit4)

fit5 = glm(recycled.objects ~ ., data = dist.data, family = "poisson")
summary(fit5)
pR2(fit5)

anova(fit1, fit2, fit3, fit4, fit5, test = "LR")

tapply(spaggs.rpj$slope_deg, spaggs.rpj$recycled.objects, summary)
tapply(spaggs.rpj$slope_deg, spaggs.rpj$recycled.objects, IQR)
pairwise.wilcox.test(spaggs.rpj$slope_deg, spaggs.rpj$recycled.objects, p.adjust.method = "bonf")
pairwise.wilcox.test(spaggs.rpj$slope_deg, spaggs.rpj$Count, p.adjust.method = "bonf")

#TODO compare with weathering stages, aggregate weight, tool numbers, etc.
comp.data = spaggs.rpj@data %>% 
        dplyr::select(recycled.objects,
                      Number_strongly_weathered_objects, 
                      Number_mildly_weathered_objects,
                      Number_weakly_weathered_objects,
                      Weight, 
                      Tool_count...40, 
                      Core_count...41, 
                      Count, 
                      dem, 
                      slope_deg, 
                      dist_to_p1, 
                      dist_to_p2, 
                      dist_to_p5) %>%
  rename(tool_count = Tool_count...40, 
         core_count = Core_count...41, 
         object_count = Count)


comp.data[is.na(comp.data)] = 0
comp.data$dist_to_sa = apply(comp.data[,c(11:13)], 1, FUN = min)

pairwise.wilcox.test(comp.data$Number_strongly_weathered_objects, comp.data$recycled.objects, p.adjust.method = "bonf")
pairwise.wilcox.test(comp.data$Number_mildly_weathered_objects, comp.data$recycled.objects, p.adjust.method = "bonf")
pairwise.wilcox.test(comp.data$Number_weakly_weathered_objects, comp.data$recycled.objects, p.adjust.method = "bonf")
#pairwise.wilcox.test(comp.data$Number_unweathered_objects, comp.data$recycled.objects, p.adjust.method = "bonf")

fit6 = glm(recycled.objects ~ Number_strongly_weathered_objects + Number_mildly_weathered_objects +
             Number_weakly_weathered_objects, data = comp.data, family = "poisson")
summary(fit6)

fit7 = glm(recycled.objects ~ object_count + tool_count + core_count, data = comp.data, family = "poisson")
summary(fit7)

fit8 = glm(recycled.objects ~ Weight + object_count + tool_count + core_count +
             Number_strongly_weathered_objects + Number_mildly_weathered_objects +
             Number_weakly_weathered_objects, data = comp.data, family = "poisson")
summary(fit8)
pR2(fit8)

fit9 = glm(recycled.objects ~ object_count + tool_count + core_count +
             Number_strongly_weathered_objects + Number_mildly_weathered_objects +
             Number_weakly_weathered_objects, data = comp.data, family = "poisson")

anova(fit6, fit7, fit8, fit9, test = "LR")

fit10 = glm(recycled.objects ~ ., data = comp.data, family = "poisson")
summary(fit10)
fit11 = glm(recycled.objects ~ dem + slope_deg + dist_to_sa + dist_to_p5 + dist_to_p2 + dist_to_p1
              + object_count + tool_count + core_count +
              Number_strongly_weathered_objects + Number_mildly_weathered_objects +
              Number_weakly_weathered_objects, data = comp.data, family = "poisson")


fit12 = glm(recycled.objects ~ dem + slope_deg + dist_to_sa 
            + object_count + tool_count + core_count +
              Number_strongly_weathered_objects + Number_mildly_weathered_objects +
              Number_weakly_weathered_objects, data = comp.data, family = "poisson")

anova(fit5, fit10, fit11, fit12, test = "LR")

summary(fit10)
pR2(fit10)

summary(fit12)
pR2(fit12)
