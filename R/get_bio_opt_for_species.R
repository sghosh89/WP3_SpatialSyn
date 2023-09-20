#rm(list=ls())
library(tidyverse)
library(dplyr)
library(here)
library(gridExtra)
library(raster)
`%notin%` <- Negate(`%in%`)
#---------- we want to get the bioopt for those species with resolved phylogeny----
df<-read.csv(here("DATA/BirdTree/species_0_250km_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)

dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_with_speciestraits_mass.csv"))
dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
df<-left_join(df,dft,by="ScientificName")

# remove the duplicated entries from df$new_BT column
df<-df%>%distinct(newBT,.keep_all = T)#253 sp.

dfm<-read.csv(here("DATA/birds_gbif_data/cleaned/birds_occurrence_metadata_added_manual.csv"))
dfm<-dfm%>%distinct(ScientificName,.keep_all = T)
dfm<-dfm%>%dplyr::select(-BirdTreeName)
df<-left_join(df,dfm,by=c("ScientificName"))


length(unique(df$ScientificName))# 253 unique species
length(unique(df$BirdTreeName))# 253 unique species

st<-unique(df$ScientificName)
bt<-unique(df$BirdTreeName)

id<-which(df$ScientificName!=df$BirdTreeName)
dfn<-df[id,]

# check if you have cleaned records for species with 
# df$ScientificName  or df$BirdTreeName

df$is_clean_csv<-NA

for(i in 1:nrow(df)){
  st<-df$ScientificName[i]
  resloc<-here("DATA/birds_gbif_data/cleaned")
  getfile<-file.exists(paste(resloc,"/",st,".csv",sep=""))
  bt<-df$BirdTreeName[i]
  if(getfile==T){
    df$is_clean_csv[i]<-1 # read from ScientificName
  }else if(file.exists(paste(resloc,"/",bt,".csv",sep=""))){
    df$is_clean_csv[i]<-2 # read from BirdTree
  }else{
    df$is_clean_csv[i]<-3 # you need to get the file for this species
  }
}
# ok, so we have all species files there: 4 sp to be read from BirdTreeName
# they are: AOU = 4120, 4881, 5275, 5660

#============get bio (1-12, 9-10, 17-18) variables for those species=====================================
bio1_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio1_1981-2010_V.2.1.tif")
rst_bio1 <- raster::raster(bio1_path)
bio12_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio12_1981-2010_V.2.1.tif")
rst_bio12 <- raster::raster(bio12_path)

#bio9_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio9_1981-2010_V.2.1.tif")
#rst_bio9 <- raster::raster(bio9_path)
#bio10_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio10_1981-2010_V.2.1.tif")
#rst_bio10 <- raster::raster(bio10_path)

#bio17_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio17_1981-2010_V.2.1.tif")
#rst_bio17 <- raster::raster(bio17_path)
#bio18_path<- here("DATA/CHELSA_v2/climatologies/1981-2010/bio/CHELSA_bio18_1981-2010_V.2.1.tif")
#rst_bio18 <- raster::raster(bio18_path)

# annual
df$bio1<-NA # See CHELSA description
df$bio12<-NA

#QUARTER
#df$bio9<-NA
#df$bio10<-NA

#df$bio17<-NA
#df$bio18<-NA


for(i in 1:nrow(df)){
  
  
  if(df$is_clean_csv[i]==1){
    sp<-df$ScientificName[i]
  }else{
    sp<-df$BirdTreeName[i]
  }
  
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
  
  print(i)
}

# scale them
df$bio1<-(df$bio1*0.1)-273.15
df$bio12<-(df$bio12*0.1)
df<-df%>%dplyr::select(-is_clean_csv)
write.csv(df,here("RESULTS/df_abund_climate_spatsyn_0_250km_with_optimal_biovar.csv"),row.names = F)

#----------------------------------------------------




