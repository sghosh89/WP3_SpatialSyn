###############################################################################################################
# This file has to be in the same folder where you save the environmental data 
# that you will extract
###############################################################################################################
rm(list=ls())
path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
library(raster)
library(sp)
rastlist <- list.files(path = "./", pattern='.tif$', all.files=TRUE, full.names=FALSE)
id<-as.integer(substr(rastlist,14,17))
idsel<-which(id%in%c(1979:2019))
rastlist<-rastlist[idsel]

allrasters <- raster::stack(rastlist)
df_lonlat<-read.csv("../uRID_lonlat_stratum.csv")

df_lonlat_table<-df_lonlat

# now make it as sp object
coordinates(df_lonlat) <- ~Longitude + Latitude
proj4string(df_lonlat) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
class(df_lonlat)

get_values<- raster::extract(allrasters, df_lonlat)
df_env<- as.data.frame(get_values)
cnm<-colnames(df_env)
cnm<-substr(cnm,8,17)
colnames(df_env)<-cnm # month_year
df_env <- cbind(df_lonlat_table,df_env)
write.csv(df_env,"../pr_monthlyvalues_extracted_for_uRID_WP3.csv",row.names=F)


