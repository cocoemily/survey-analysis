## DATA CLEANING
library(tidyverse)
library(rgdal)
library(raster)
library(sf)
library(sfhotspot)
library(spDataLarge)
library(spdep)
library(tmap)
library(spatstat)

source("bordian_types_dictionary.R")

#### CREATE NUMERICS ####
data$Artfct_t = factor(data$Artfct_t, levels = c("complete_flake", "broken_flake", "tool", "tool_fragment", "core", "core_fragment", "shatter"))
data$atype.num = unclass(data$Artfct_t)

data$Wthrng_c = factor(data$Wthrng_c, levels = c("strongly_weathered", "mildly_weathered", "weakly_weathered", "not_weathered", "other"))
data$wt.num = unclass(data$Wthrng_c)

data$flake.type = ifelse(str_detect(data$Flk_ty, pattern = "flake"), "flake", data$Flk_ty)
data$flake.type = factor(data$flake.type, levels = c("flake", "blade", "bladelet", "other"))
data$ftype.num = unclass(data$flake.type)

data$Flk_tr = factor(data$Flk_tr, levels = c("feather", "hinge", "plunge", "step", "other"))
data$flk.tr.num = unclass(data$Flk_tr)

data$retouch.side = ifelse(str_detect(data$Rtch_s, pattern = " "), "bifacial", data$Rtch_s)
data$retouch.side = factor(data$retouch.side, levels = c("dorsal", "ventral", "bifacial"))
data$rets.num = unclass(data$retouch.side)

#data$plat.type 
#need to think about how to do this one

data$tool.type = ifelse(str_detect(data$Tl_ty, "notch denticulate"), "notch/denticulate", 
                        ifelse(str_detect(data$Tl_ty, " "), "multiple", 
                               data$Tl_ty))
data$tool.type = factor(data$tool.type, levels = c("notch", "denticulate", "notch/denticulate", "scraper", "point", "biface", "multiple", "other"))
data$ttype.num = unclass(data$tool.type)

data$Bordian_name = ""
for(i in 1:nrow(data)) {
  #print(data$Brdn_[i])
  if(!is.na(data$Brdn_[i])) {
    data$Bordian_name[i] = bordian_types[data$Brdn_[i]]
  }
}
data$Bordian_name = factor(data$Bordian_name, levels = bordian_levels)
data$bord.num = unclass(data$Bordian_name)



#### list of numerics ####
# 1) data$Pltfrm_w
# 2) data$Pltfrm_th
# 3) data$Drsl_flk_s
# 4) data$Retch ##true/false
# 5) data$Edg_d ##true/false
# 6) data$Crtx_p
# 7) data$Weght
# 8) data$Lngth
# 9) data$Width
# 10) data$Thckn
# 11) data$rcycl ##true/false
# 12) data$dbl_p ##true/false
# 13) data$atype.num
# 14) data$wt.num
# 15) data$ftype.num
# 16) data$flk.tr.num
# 17) data$rets.num
# 18) data$ttype.num
# 19) data$bord.num

#presence data
data$compl_flk = ifelse(data$Artfct_t == "complete_flake", 1, 0)
data$broke_flk = ifelse(data$Artfct_t == "broken_flake", 1, 0)
data$tool = ifelse(data$Artfct_t == "tool", 1, 0)
data$tool_frag = ifelse(data$Artfct_t == "tool_fragment", 1, 0)
data$core = ifelse(data$Artfct_t == "core", 1, 0)
data$core_frag = ifelse(data$Artfct_t == "core_fragment", 1, 0)
data$shatter = ifelse(data$Artfct_t == "shatter", 1, 0)

data$str_weather = ifelse(data$Wthrng_c == "strongly_weathered", 1, 0)
data$mid_weather = ifelse(data$Wthrng_c == "mildly_weathered", 1, 0)
data$weak_weather = ifelse(data$Wthrng_c == "weakly_weathered", 1, 0)
data$not_weather = ifelse(data$Wthrng_c == "not_weathered", 1, 0)
data$oth_weather = ifelse(data$Wthrng_c == "other", 1, 0)

data$type_flake = ifelse(data$flake.type == "flake", 1, 0)
data$type_blade = ifelse(data$flake.type == "blade", 1, 0)
data$type_bladelet = ifelse(data$flake.type == "bladelet", 1, 0)
data$type_oth = ifelse(data$flake.type == "other", 1, 0)

data$rside_dorsal = ifelse(data$retouch.side == "dorsal", 1, 0)
data$rside_ventral = ifelse(data$retouch.side == "ventral", 1, 0)
data$rside_bifacial = ifelse(data$retouch.side == "bifacial", 1, 0)

data$ttype_notch = ifelse(data$tool.type == "notch", 1, 0)
data$ttype_dent = ifelse(data$tool.type == "denticulate", 1, 0)
data$ttype_nd = ifelse(data$tool.type == "notch/denticulate", 1, 0)
data$ttype_scraper = ifelse(data$tool.type == "scraper", 1, 0)
data$ttype_point = ifelse(data$tool.type == "point", 1, 0)
data$ttype_biface = ifelse(data$tool.type == "biface", 1, 0)
data$ttype_mult = ifelse(data$tool.type == "multiple", 1, 0)
data$ttype_oth = ifelse(data$tool.type == "other", 1, 0)