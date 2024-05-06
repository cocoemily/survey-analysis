set.seed(120109)

library(tidyverse)
library(rgdal)
library(raster)
library(sf)
library(spatstat)
library(rasterVis)
library(RColorBrewer)
library(ggpubr)

p2.artifacts = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/data/artifact-shapefiles/", layer = "p2-artifacts")
p2.window = readOGR("/Users/emilycoco/Desktop/NYU/Dissertation-Research/Survey/survey-analysis/output/square-extents/", layer = "p2-sqs")

data = p2.artifacts
window = p2.window

source("intrasite-analysis-scripts/clean-data.R")
colpal = rev(c("#E69F00", "#0C7BDC", "#009E73", "grey80"))

st_data = st_transform(st_as_sf(data), 32642) #WGS 84 / UTM zone 42N
win = as.owin(st_transform(st_as_sf(window), 32642))


artifact.ppp = as.ppp(st_data)
marks(artifact.ppp) = NULL
Window(artifact.ppp) = win

rcycl.ppp = as.ppp(st_data %>% filter(rcycl == 1))
marks(rcycl.ppp) = NULL
Window(rcycl.ppp) = win

nrcycl.ppp = as.ppp(st_data %>% filter(rcycl == 0))
marks(nrcycl.ppp) = NULL
Window(nrcycl.ppp) = win


LR = scanLRTS(rcycl.ppp, r = 2*bw.diggle(artifact.ppp))
pvals = eval.im(pchisq(LR, df = 1, lower.tail=F))
rcycl.hs = as.im(pvals < 0.05)
Window(rcycl.hs) = win
plot(rcycl.hs)

nLR = scanLRTS(nrcycl.ppp, r = 2*bw.diggle(artifact.ppp))
npvals = eval.im(pchisq(nLR, df = 1, lower.tail=F))
nrcycl.hs = as.im(npvals < 0.05)
Window(nrcycl.hs) = win
plot(nrcycl.hs)

#### Artifact types Hotspots Exploration ####
aLR = scanLRTS(artifact.ppp, r = 2*bw.diggle(artifact.ppp))
aPvals = eval.im(pchisq(aLR, df = 1, lower.tail=F))
artifact.hs = as.im(aPvals < 0.05)
Window(artifact.hs) = win
plot(artifact.hs)

# cfLR = scanLRTS(compl_flk.ppp, r = 2*bw.diggle(artifact.ppp))
# cfPvals = eval.im(pchisq(cfLR, df = 1, lower.tail=F))
# cf.hs = as.im(cfPvals < 0.05)
# Window(cf.hs) = win
# plot(cf.hs)
# 
# bfLR = scanLRTS(broke_flk.ppp, r = 2*bw.diggle(artifact.ppp))
# bfPvals = eval.im(pchisq(bfLR, df = 1, lower.tail=F))
# bf.hs = as.im(bfPvals < 0.05)
# Window(bf.hs) = win
# plot(bf.hs)

pppdata = st_data %>% filter(compl_flk == 1 | broke_flk == 1) %>% filter(rcycl == 0)
allflake.ppp = as.ppp(pppdata)
marks(allflake.ppp) = NULL
Window(allflake.ppp) = win
flakeLR = scanLRTS(allflake.ppp, r = 2*bw.diggle(artifact.ppp))
flakePvals = eval.im(pchisq(flakeLR, df = 1, lower.tail=F))
flake.hs = as.im(flakePvals < 0.05)
Window(flake.hs) = win
plot(flake.hs)

pppdata = st_data %>% filter(core == 1 | core_frag == 1) %>% filter(rcycl == 0)
allcore.ppp = as.ppp(pppdata)
marks(allcore.ppp) = NULL
Window(allcore.ppp) = win
coreLR = scanLRTS(allcore.ppp, r = 2*bw.diggle(artifact.ppp))
corePvals = eval.im(pchisq(coreLR, df = 1, lower.tail=F))
core.hs = as.im(corePvals < 0.05)
Window(core.hs) = win
plot(core.hs)

pppdata = st_data %>% filter(tool == 1 | tool_frag == 1) %>% filter(rcycl == 0)
alltool.ppp = as.ppp(pppdata)
marks(alltool.ppp) = NULL
Window(alltool.ppp) = win
toolLR = scanLRTS(alltool.ppp, r = 2*bw.diggle(artifact.ppp))
toolPvals = eval.im(pchisq(toolLR, df = 1, lower.tail=F))
tool.hs = as.im(toolPvals < 0.05)
Window(tool.hs) = win
plot(tool.hs)


#####plotting#####
colpal = rev(c("#E69F00", "#0C7BDC", "#009E73", "grey80"))

#at.hs.overlap = (artifact.hs*2) - rcycl.hs
at.hs.overlap = (nrcycl.hs*2) - rcycl.hs
#plot.im(at.hs.overlap, useRaster=F, ribn = 4)
at.hs.rast = raster(at.hs.overlap)
at.hs.rast = as.factor(at.hs.rast)
at.hs <- levels(at.hs.rast)[[1]]
at.hs[['hotspot']] = rev(c("artifacts", "overlap", "none", "recycled artifacts"))
at.hs[['order']] = (c(2,4,3,1))
at.hs = at.hs[order(-at.hs$order),]
levels(at.hs.rast) = at.hs
perc.at = as.data.frame(freq(at.hs.rast))
perc.at[['hotspot']] = c(rev(c("artifacts", "overlap", "none", "recycled artifacts")), "NA")
perc.at$percentage = perc.at$count/cellStats(at.hs.rast, function(i, ...) sum(!is.na(i)))
perc.at$recycle.percentage = perc.at$count/(perc.at$count[1] + perc.at$count[3])
perc.at$artifacts.percentage = perc.at$count/(perc.at$count[4] + perc.at$count[3])
perc.at$overlap = "artifacts"
percentage = max(round(perc.at[which(perc.at$value == 1),]$artifacts.percentage*100, 2),0)
plot1 = levelplot(at.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of non-recycled artifact hotspots"))
plot(plot1)

f.hs.overlap = (flake.hs*2) - rcycl.hs
#plot(f.hs.overlap, useRaster=F)
f.hs.rast = raster(f.hs.overlap)
f.hs.rast = as.factor(f.hs.rast)
f.hs <- levels(f.hs.rast)[[1]]
f.hs[['hotspot']] = rev(c("flakes", "overlap", "none", "recycled artifacts"))
f.hs[['order']] = (c(2,4,3,1))
f.hs = f.hs[order(-f.hs$order),]
levels(f.hs.rast) = f.hs
perc.f = as.data.frame(freq(f.hs.rast))
perc.f[['hotspot']] = c(rev(c("flakes", "overlap", "none", "recycled artifacts")), "NA")
perc.f$percentage = perc.f$count/cellStats(f.hs.rast, function(i, ...) sum(!is.na(i)))
perc.f$recycle.percentage = perc.f$count/(perc.f$count[1] + perc.f$count[3])
perc.f$artifacts.percentage = perc.f$count/(perc.f$count[4] + perc.f$count[3])
perc.f$overlap = "flakes"
percentage = max(round(perc.f[which(perc.f$value == 1),]$artifacts.percentage*100, 2),0)
plot2 = levelplot(f.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of non-recycled flake hotspots"))
plot(plot2)

# cf.hs.overlap = (cf.hs*2) - rcycl.hs
# #plot(cf.hs.overlap, useRaster=F)
# cf.hs.rast = raster(cf.hs.overlap)
# cf.hs.rast = as.factor(cf.hs.rast)
# cf.hs <- levels(cf.hs.rast)[[1]]
# cf.hs[['hotspot']] = rev(c("complete flakes", "overlap", "none", "recycled artifacts"))
# cf.hs[['order']] = (c(2,4,3,1))
# cf.hs = cf.hs[order(-cf.hs$order),]
# levels(cf.hs.rast) = cf.hs
# perc.cf = as.data.frame(freq(cf.hs.rast))
# perc.cf[['hotspot']] = c(rev(c("complete flakes", "overlap", "none", "recycled artifacts")), "NA")
# perc.cf$percentage = perc.cf$count/cellStats(cf.hs.rast, function(i, ...) sum(!is.na(i)))
# perc.cf$recycle.percentage = perc.cf$count/(perc.cf$count[1] + perc.cf$count[3])
# perc.cf$artifacts.percentage = perc.cf$count/(perc.cf$count[4] + perc.cf$count[3])
# perc.cf$overlap = "complete flakes"
# percentage = max(round(perc.cf[which(perc.cf$value == 1),]$artifacts.percentage*100, 2),0)
# plot2 = levelplot(cf.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of complete flake hotspots"))
# 
# bf.hs.overlap = (bf.hs*2) - rcycl.hs
# #plot(bf.hs.overlap, useRaster=F)
# bf.hs.rast = raster(bf.hs.overlap)
# bf.hs.rast = as.factor(bf.hs.rast)
# bf.hs <- levels(bf.hs.rast)[[1]]
# bf.hs[['hotspot']] = rev(c("broken flakes", "overlap", "none", "recycled artifacts"))
# bf.hs[['order']] = (c(2,4,3,1))
# bf.hs = bf.hs[order(-bf.hs$order),]
# levels(bf.hs.rast) = bf.hs
# perc.bf = as.data.frame(freq(bf.hs.rast))
# perc.bf[['hotspot']] = c(rev(c("broken flakes", "overlap", "none", "recycled artifacts")), "NA")
# perc.bf$percentage = perc.bf$count/cellStats(bf.hs.rast, function(i, ...) sum(!is.na(i)))
# perc.bf$recycle.percentage = perc.bf$count/(perc.bf$count[1] + perc.bf$count[3])
# perc.bf$artifacts.percentage = perc.bf$count/(perc.bf$count[4] + perc.bf$count[3])
# perc.bf$overlap = "broken flakes"
# percentage = max(round(perc.bf[which(perc.bf$value == 1),]$artifacts.percentage*100, 2), 0)
# plot3 = levelplot(bf.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of broken flake hotspots"))

# c.hs.overlap = (core.hs*2) - rcycl.hs
# plot(c.hs.overlap, useRaster=F)
# c.hs.rast = raster(c.hs.overlap)
# c.hs.rast = as.factor(c.hs.rast)
# c.hs <- levels(c.hs.rast)[[1]]
# c.hs[['hotspot']] = rev(c("cores", "overlap", "none", "recycled artifacts"))
# c.hs[['order']] = (c(2,4,3,1))
# c.hs = c.hs[order(-c.hs$order),]
# levels(c.hs.rast) = c.hs
# perc.c = as.data.frame(freq(c.hs.rast))
# perc.c[['hotspot']] = c(rev(c("cores", "overlap", "none", "recycled artifacts")), "NA")
# perc.c$percentage = perc.c$count/cellStats(c.hs.rast, function(i, ...) sum(!is.na(i)))
# perc.c$recycle.percentage = perc.c$count/(perc.c$count[1] + perc.c$count[3])
# perc.c$artifacts.percentage = perc.c$count/(perc.c$count[4] + perc.c$count[3])
# perc.c$overlap = "cores"
# percentage = max(round(perc.c[which(perc.c$value == 1),]$artifacts.percentage*100, 2),0)
# plot4 = levelplot(c.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of core hotspots"))

t.hs.overlap = (tool.hs*2) - rcycl.hs
#plot(t.hs.overlap, useRaster=F)
t.hs.rast = raster(t.hs.overlap)
t.hs.rast = as.factor(t.hs.rast)
t.hs <- levels(t.hs.rast)[[1]]
t.hs[['hotspot']] = rev(c("tools", "overlap", "none", "recycled artifacts"))
t.hs[['order']] = (c(2,4,3,1))
t.hs = t.hs[order(-t.hs$order),]
levels(t.hs.rast) = t.hs
perc.t = as.data.frame(freq(t.hs.rast))
perc.t[['hotspot']] = c(rev(c("tools", "overlap", "none", "recycled artifacts")), "NA")
perc.t$percentage = perc.t$count/cellStats(t.hs.rast, function(i, ...) sum(!is.na(i)))
perc.t$recycle.percentage = perc.t$count/(perc.t$count[1] + perc.t$count[3])
perc.t$artifacts.percentage = perc.t$count/(perc.t$count[4] + perc.t$count[3])
perc.t$overlap = "tools"
percentage = max(round(perc.t[which(perc.t$value == 1),]$artifacts.percentage*100, 2),0)
plot5 = levelplot(t.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of non-recycled tool hotspots"))
plot(plot5)

allperc.a = rbind(perc.at, perc.f, perc.t) %>% filter(value == 1)

tiff(filename = "figures/p2/artifact-recycling-overlap.tiff", 
     units="in", 
     width=5*2, 
     height=4*2, 
     pointsize=12, 
     res=300)
plot1
dev.off()

tiff(filename = "figures/p2/flakes-recycling-overlap.tiff",
     units="in",
     width=5*2,
     height=4*2,
     pointsize=12,
     res=300)
plot2
dev.off()

# tiff(filename = "figures/p2/complete_flakes-recycling-overlap.tiff", 
#      units="in", 
#      width=5*2, 
#      height=4*2, 
#      pointsize=12, 
#      res=300)
# plot2
# dev.off()
# 
# tiff(filename = "figures/p2/broken_flakes-recycling-overlap.tiff", 
#      units="in", 
#      width=5*2, 
#      height=4*2, 
#      pointsize=12, 
#      res=300)
# plot3
# dev.off()
# 
# tiff(filename = "figures/p2/cores-recycling-overlap.tiff", 
#      units="in", 
#      width=5*2, 
#      height=4*2, 
#      pointsize=12, 
#      res=300)
# plot4
# dev.off()

tiff(filename = "figures/p2/tools-recycling-overlap.tiff", 
     units="in", 
     width=5*2, 
     height=4*2, 
     pointsize=12, 
     res=300)
plot5
dev.off()

#### WEATHERING CLASSES Hotspots ####
pppdata = st_data %>% filter(str_weather == 1) %>% filter(rcycl == 0)
if(nrow(pppdata) > 0) {
  str_weather.ppp = as.ppp(pppdata)
  marks(str_weather.ppp) = NULL
  Window(str_weather.ppp) = win
  str_weather.dens = density(str_weather.ppp, sigma = bw.diggle, adjust = 2)
}
swLR = scanLRTS(str_weather.ppp, r = 2*bw.diggle(artifact.ppp))
pvals = eval.im(pchisq(swLR, df = 1, lower.tail=F))
strw.hs = as.im(pvals < 0.05)
Window(strw.hs) = win
plot(strw.hs)


pppdata = st_data %>% filter(mid_weather == 1) %>% filter(rcycl == 0)
if(nrow(pppdata) > 0) {
  mid_weather.ppp = as.ppp(pppdata)
  marks(mid_weather.ppp) = NULL
  Window(mid_weather.ppp) = win
  mid_weather.dens = density(mid_weather.ppp, sigma = bw.diggle, adjust = 2)
}
mwLR = scanLRTS(mid_weather.ppp, r = 2*bw.diggle(artifact.ppp))
pvals = eval.im(pchisq(mwLR, df = 1, lower.tail=F))
midw.hs = as.im(pvals < 0.05)
Window(midw.hs) = win
plot(midw.hs)


pppdata = st_data %>% filter(weak_weather == 1) %>% filter(rcycl == 0)
if(nrow(pppdata) > 0) {
  weak_weather.ppp = as.ppp(pppdata)
  marks(weak_weather.ppp) = NULL
  Window(weak_weather.ppp) = win
  weak_weather.dens = density(weak_weather.ppp, sigma = bw.diggle, adjust = 2)
}
wwLR = scanLRTS(weak_weather.ppp, r = 2*bw.diggle(artifact.ppp))
pvals = eval.im(pchisq(wwLR, df = 1, lower.tail=F))
weakw.hs = as.im(pvals < 0.05)
Window(weakw.hs) = win
plot(weakw.hs)

pppdata = st_data %>% filter(not_weather == 1) %>% filter(rcycl == 0)
if(nrow(pppdata) > 0) {
  not_weather.ppp = as.ppp(pppdata)
  marks(not_weather.ppp) = NULL
  Window(not_weather.ppp) = win
  not_weather.dens = density(not_weather.ppp, sigma = bw.diggle, adjust = 2)
}
nwLR = scanLRTS(not_weather.ppp, r = 2*bw.diggle(artifact.ppp))
pvals = eval.im(pchisq(nwLR, df = 1, lower.tail=F))
notw.hs = as.im(pvals < 0.05)
Window(notw.hs) = win
plot(notw.hs)

#####plotting#####
sw.hs.overlap = (strw.hs*2) - rcycl.hs
plot.im(sw.hs.overlap, useRaster=F)
sw.hs.rast = raster(sw.hs.overlap)
sw.hs.rast = as.factor(sw.hs.rast)
sw.hs <- levels(sw.hs.rast)[[1]]
sw.hs[['hotspot']] = rev(c("strongly weathered", "overlap", "none", "recycled artifacts"))
sw.hs[['order']] = (c(2,4,3,1))
sw.hs = sw.hs[order(-sw.hs$order),]
levels(sw.hs.rast) = sw.hs
perc.sw = as.data.frame(freq(sw.hs.rast))
perc.sw[['hotspot']] = c(rev(c("strongly weathered", "overlap", "none", "recycled artifacts")), "NA")
perc.sw$percentage = perc.sw$count/cellStats(sw.hs.rast, function(i, ...) sum(!is.na(i)))
perc.sw$recycle.percentage = perc.sw$count/(perc.sw$count[1] + perc.sw$count[3])
perc.sw$artifacts.percentage = perc.sw$count/(perc.sw$count[4] + perc.sw$count[3])
perc.sw$overlap = "strongly weathered"
percentage = max(round(perc.sw[which(perc.sw$value == 1),]$artifacts.percentage*100, 2),0)
plot1 = levelplot(sw.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of strongly weathered NR artifact hotspots"))
plot(plot1)

mw.hs.overlap = (midw.hs*2) - rcycl.hs
plot.im(mw.hs.overlap, useRaster=F)
mw.hs.rast = raster(mw.hs.overlap)
mw.hs.rast = as.factor(mw.hs.rast)
mw.hs <- levels(mw.hs.rast)[[1]]
mw.hs[['hotspot']] = rev(c("mildly weathered", "overlap", "none", "recycled artifacts"))
mw.hs[['order']] = (c(2,4,3,1))
mw.hs = mw.hs[order(-mw.hs$order),]
levels(mw.hs.rast) = mw.hs
perc.mw = as.data.frame(freq(mw.hs.rast))
perc.mw[['hotspot']] = c(rev(c("mildly weathered", "overlap", "none", "recycled artifacts")), "NA")
perc.mw$percentage = perc.mw$count/cellStats(mw.hs.rast, function(i, ...) sum(!is.na(i)))
perc.mw$recycle.percentage = perc.mw$count/(perc.mw$count[1] + perc.mw$count[3])
perc.mw$artifacts.percentage = perc.mw$count/(perc.mw$count[4] + perc.mw$count[3])
perc.mw$overlap = "mildly weathered"
percentage = max(round(perc.mw[which(perc.mw$value == 1),]$artifacts.percentage*100, 2),0)
plot2 = levelplot(mw.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of mildly weathered NR artifact hotspots"))
plot(plot2)

ww.hs.overlap = (weakw.hs*2) - rcycl.hs
plot.im(ww.hs.overlap, useRaster=F)
ww.hs.rast = raster(ww.hs.overlap)
ww.hs.rast = as.factor(ww.hs.rast)
ww.hs <- levels(ww.hs.rast)[[1]]
ww.hs[['hotspot']] = rev(c("weakly weathered", "overlap", "none", "recycled artifacts"))
ww.hs[['order']] = (c(2,4,3,1))
ww.hs = ww.hs[order(-ww.hs$order),]
levels(ww.hs.rast) = ww.hs
perc.ww = as.data.frame(freq(ww.hs.rast))
perc.ww[['hotspot']] = c(rev(c("weakly weathered", "overlap", "none", "recycled artifacts")), "NA")
perc.ww$percentage = perc.ww$count/cellStats(ww.hs.rast, function(i, ...) sum(!is.na(i)))
perc.ww$recycle.percentage = perc.ww$count/(perc.ww$count[1] + perc.ww$count[3])
perc.ww$artifacts.percentage = perc.ww$count/(perc.ww$count[4] + perc.ww$count[3])
perc.ww$overlap = "weakly weathered"
percentage = max(round(perc.ww[which(perc.ww$value == 1),]$artifacts.percentage*100, 2),0)
plot3 = levelplot(ww.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of weakly weathered NR artifact hotspots"))
plot(plot3)

nw.hs.overlap = (notw.hs*2) - rcycl.hs
plot.im(nw.hs.overlap, useRaster=F)
nw.hs.rast = raster(nw.hs.overlap)
nw.hs.rast = as.factor(nw.hs.rast)
nw.hs <- levels(nw.hs.rast)[[1]]
nw.hs[['hotspot']] = rev(c("unweathered", "overlap", "none", "recycled artifacts"))
nw.hs[['order']] = (c(2,4,3,1))
nw.hs = nw.hs[order(-nw.hs$order),]
levels(nw.hs.rast) = nw.hs
perc.nw = as.data.frame(freq(nw.hs.rast))
perc.nw[['hotspot']] = c(rev(c("unweathered","overlap", "none", "recycled artifacts")), "NA")
perc.nw$percentage = perc.nw$count/cellStats(nw.hs.rast, function(i, ...) sum(!is.na(i)))
perc.nw$recycle.percentage = perc.nw$count/(perc.nw$count[1] + perc.nw$count[3])
perc.nw$artifacts.percentage = perc.nw$count/(perc.nw$count[4] + perc.nw$count[3])
perc.nw$overlap = "unweathered"
percentage = max(round(perc.nw[which(perc.nw$value == 1),]$artifacts.percentage*100, 2),0)
plot4 = levelplot(nw.hs.rast, col.regions = colpal, main = paste0("recycled artifact hotspots overlap with ", percentage, "% of unweathered NR artifact hotspots"))
plot(plot4)

allperc = rbind(perc.sw, perc.mw, perc.ww, perc.nw) %>% filter(value == 1)

tiff(filename = "figures/p2/strongly-weathered-recycling-overlap.tiff", 
     units="in", 
     width=5*2, 
     height=4*2, 
     pointsize=12, 
     res=300)
plot1
dev.off()

tiff(filename = "figures/p2/mildly-weathered-recycling-overlap.tiff", 
     units="in", 
     width=5*2, 
     height=4*2, 
     pointsize=12, 
     res=300)
plot2
dev.off()

tiff(filename = "figures/p2/weakly-weathered-recycling-overlap.tiff", 
     units="in", 
     width=5*2, 
     height=4*2, 
     pointsize=12, 
     res=300)
plot3
dev.off()

tiff(filename = "figures/p2/unweathered-recycling-overlap.tiff", 
     units="in", 
     width=5*2, 
     height=4*2, 
     pointsize=12, 
     res=300)
plot4
dev.off()
