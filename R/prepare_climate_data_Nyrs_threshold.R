# first get climate data in required format
# we will first save the climate data for sites that were selected for each species spat syn computation
# this is more flexible function than before: included yr_threshold argument

#rm(list=ls())
library(here)
library(dplyr)
library(tidyverse)
library(maps)
library(pracma)

prepare_climate_data_Nyrs_threshold<-function(climvar="tas", yr_threshold=40){
  
  cdat<-read.csv(here(paste("DATA/CHELSA_v2/monthly/",climvar,"_monthlyvalues_extracted_for_uRID_WP3.csv",sep="")))
  cdat<-cdat%>%dplyr::select(-Longitude, -Latitude, -Stratum)
  longcdat<- gather(cdat,key="key",value="value", -uRID)
  longcdat<-longcdat%>%separate(key,c("var","month","year"))
  longcdat$year<-as.integer(longcdat$year)
  #longcdat<-longcdat%>%filter(year%in%c(1997:2019))# remove year=1996 
  
  longcdat<-longcdat%>%group_by(uRID,year)%>%summarise(meanval=mean(value,na.rm=T))%>%ungroup()
  zz<-longcdat%>%spread(year, meanval)
  cdat<-t(zz)
  colnames(cdat)<-cdat[1,]
  cdat<-cdat[-1,]
  dim(cdat)# year by site
  class(cdat[1,])# character
  rnm<-rownames(cdat)
  cnm<-colnames(cdat)
  cdat<-matrix(as.numeric(cdat),    # Convert from character matrix to numeric matrix
               ncol = ncol(cdat))
  rownames(cdat)<-rnm
  colnames(cdat)<-cnm
  if(climvar=="pr"){
    cdat<-rbind(cdat,NA*numeric(ncol(cdat)))
    rownames(cdat)[nrow(cdat)]<-"2019"
  }
  saveRDS(cdat,here(paste("RESULTS/year_by_site_",climvar,".RDS",sep="")))
  
  
  # These species are chosen
  if(yr_threshold==40){
    df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))
  }else{
    df<-read.csv(here(paste("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites_min",yr_threshold,"yr.csv",sep="")))
  }
  
  for(j in 1:nrow(df)){
    #j<-1
    givenAOU<-df$AOU[j]
    
    if(yr_threshold==40){
      resloc<-here(paste("RESULTS/AOU_",givenAOU,sep=""))
    }else{
      resloc<-here(paste("RESULTS/yr_threshold_",yr_threshold,"/AOU_",givenAOU,sep=""))
    }
    
    distm_sel<-readRDS(here(paste(resloc,"/distm_sel.RDS",sep="")))
    
    selsite<-rownames(distm_sel)
    id<-which(colnames(cdat)%in%selsite)
    cdat_sel<-cdat[,id] # climate timeseries selected for chosen sites
    
    mat<-as.matrix(cdat_sel)
    dmat<-pracma::detrend(cdat_sel) # detrended climate timeseries, removing linear trend from each column
    
    d_allsite<-vector(mode="list",length=ncol(mat))
    names(d_allsite)<-colnames(mat)
    d_allsite_detrend<-d_allsite
    
    for(i in 1:length(d_allsite)){
      d_allsite[[i]]<-data.frame(Year=rownames(mat),Dat=mat[,i])
      d_allsite_detrend[[i]]<-data.frame(Year=rownames(dmat),Dat=dmat[,i])
    }
    
    saveRDS(mat,paste(resloc,"/year_by_site_",climvar,".RDS",sep=""))
    saveRDS(dmat,paste(resloc,"/year_by_site_",climvar,"_detrended.RDS",sep=""))
    
    saveRDS(d_allsite,paste(resloc,"/",climvar,"_data_selectedsitelist.RDS",sep=""))
    saveRDS(d_allsite_detrend,paste(resloc,"/",climvar,"_detrended_data_selectedsitelist.RDS",sep=""))
    print(j)
  }
  cat("---- done ----\n")
}

#call the function
climvar<-"pr"
prepare_climate_data_Nyrs_threshold(climvar=climvar,yr_threshold = 32)
prepare_climate_data_Nyrs_threshold(climvar=climvar,yr_threshold = 36)
prepare_climate_data_Nyrs_threshold(climvar=climvar,yr_threshold = 40)

climvar<-"tas"
prepare_climate_data_Nyrs_threshold(climvar=climvar,yr_threshold = 32)
prepare_climate_data_Nyrs_threshold(climvar=climvar,yr_threshold = 36)
prepare_climate_data_Nyrs_threshold(climvar=climvar,yr_threshold = 40)

climvar<-"tasmax"
prepare_climate_data_Nyrs_threshold(climvar=climvar,yr_threshold = 32)
prepare_climate_data_Nyrs_threshold(climvar=climvar,yr_threshold = 36)
prepare_climate_data_Nyrs_threshold(climvar=climvar,yr_threshold = 40)


