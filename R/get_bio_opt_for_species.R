#rm(list=ls())
library(tidyverse)
library(dplyr)
library(here)
library(gridExtra)
library(raster)
#--------------------
df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km.csv"))
splist<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.csv"))
splist<-splist%>%dplyr::select(AOU,ScientificName)
df<-left_join(df,splist,by="AOU")

dftrait<-read.csv(here("DATA/AVONET/bird_traits_from_AVONET.csv"))
dftrait<-dftrait%>%dplyr::select(spname,possible_sp)
dftrait<-dftrait%>%distinct(spname,.keep_all = T)
df<-left_join(df,dftrait,by=c("ScientificName"="spname"))

dfm<-read.csv(here("DATA/birds_gbif_data/cleaned/birds_occurrence_metadata_added_manual.csv"))
dfm<-dfm%>%distinct(ScientificName,.keep_all = T)
df<-left_join(df,dfm,by=c("ScientificName"))

# possible_sp column contains species for which you need the cleaned gbif data

length(unique(df$possible_sp))# 257 unique species

length(unique(df$ScientificName))# 262 unique species

#============get bio (1-12, 9-10, 17-18) variables for those species=====================================
bio1_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio1_1981-2010_V.2.1.tif")
rst_bio1 <- raster::raster(bio1_path)
bio12_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio12_1981-2010_V.2.1.tif")
rst_bio12 <- raster::raster(bio12_path)

bio9_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio9_1981-2010_V.2.1.tif")
rst_bio9 <- raster::raster(bio9_path)
bio10_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio10_1981-2010_V.2.1.tif")
rst_bio10 <- raster::raster(bio10_path)

bio17_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio17_1981-2010_V.2.1.tif")
rst_bio17 <- raster::raster(bio17_path)
bio18_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio18_1981-2010_V.2.1.tif")
rst_bio18 <- raster::raster(bio18_path)

# annual
df$bio1<-NA # See CHELSA description
df$bio12<-NA

#QUARTER
df$bio9<-NA
df$bio10<-NA

df$bio17<-NA
df$bio18<-NA


for(i in 1:nrow(df)){
  
  sp<-df$possible_sp[i]
 # temp<-file.exists(here(paste("DATA/birds_gbif_data/cleaned/",sp,".csv",sep="")))
  
  # read the species occurence data
  occ_dat<-read.csv(here(paste("DATA/birds_gbif_data/cleaned/",sp,".csv",sep="")))
  df_lonlat_table<-occ_dat%>%dplyr::select(decimalLongitude,decimalLatitude)
  
  # now make it as sp object
  coordinates(df_lonlat_table) <- ~decimalLongitude + decimalLatitude
  proj4string(df_lonlat_table) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  #---------- for bio1-----------------------------
  get_values_bio1<- raster::extract(rst_bio1, df_lonlat_table)
  
  occ_dat$bio1<-get_values_bio1
  
  xx<-occ_dat
  xxbio1<-xx%>%dplyr::select(species,year,bio1)%>%distinct(bio1,.keep_all = T)
  #xxbio1<-xxbio1%>%filter(year<1980)
  xxbio1<-xxbio1%>%arrange(desc(bio1))
  df$bio1[i]<-mean(xxbio1$bio1,na.rm=T)
  
  #-------------- for bio12 ---------------------
  get_values_bio12<- raster::extract(rst_bio12, df_lonlat_table)
  
  occ_dat$bio12<-get_values_bio12
  
  xx<-occ_dat
  xxbio12<-xx%>%dplyr::select(species,year,bio12)%>%distinct(bio12,.keep_all = T)
  #xxbio12<-xxbio12%>%arrange(desc(bio12))
  df$bio12[i]<-mean(xxbio12$bio12,na.rm=T)
  
  #-------------- for bio9 ---------------------
  get_values_bio9<- raster::extract(rst_bio9, df_lonlat_table)
  
  occ_dat$bio9<-get_values_bio9
  
  xx<-occ_dat
  xxbio9<-xx%>%dplyr::select(species,year,bio9)%>%distinct(bio9,.keep_all = T)
  #xxbio9<-xxbio9%>%arrange(desc(bio9))
  df$bio9[i]<-mean(xxbio9$bio9,na.rm=T)
  
  #-------------- for bio10 ---------------------
  get_values_bio10<- raster::extract(rst_bio10, df_lonlat_table)
  
  occ_dat$bio10<-get_values_bio10
  
  xx<-occ_dat
  xxbio10<-xx%>%dplyr::select(species,year,bio10)%>%distinct(bio10,.keep_all = T)
  df$bio10[i]<-mean(xxbio10$bio10,na.rm=T)
  
  #-------------- for bio17 ---------------------
  get_values_bio17<- raster::extract(rst_bio17, df_lonlat_table)
  
  occ_dat$bio17<-get_values_bio17
  
  xx<-occ_dat
  xxbio17<-xx%>%dplyr::select(species,year,bio17)%>%distinct(bio17,.keep_all = T)
  df$bio17[i]<-mean(xxbio17$bio17,na.rm=T)
  
  #-------------- for bio18 ---------------------
  get_values_bio18<- raster::extract(rst_bio18, df_lonlat_table)
  
  occ_dat$bio18<-get_values_bio18
  
  xx<-occ_dat
  xxbio18<-xx%>%dplyr::select(species,year,bio18)%>%distinct(bio18,.keep_all = T)
  df$bio18[i]<-mean(xxbio18$bio18,na.rm=T)
  
  print(i)
}

# scale them
df$bio1<-(df$bio1*0.1)-273.15
df$bio12<-(df$bio12*0.1)

df$bio9<-(df$bio9*0.1)-273.15
df$bio10<-(df$bio10*0.1)-273.15

df$bio17<-(df$bio17*0.1)
df$bio18<-(df$bio18*0.1)

write.csv(df,here("RESULTS/df_abund_climate_spatsyn_0_250km_with_optimal_biovar.csv"),row.names = F)

#----------------------------------------------------




