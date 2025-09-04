library(here)
library(tidyverse)

xweather<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/weather.csv")) # weather data
colnames(xweather)
range(xweather$Year)
unique(xweather$Month)

tab<-xweather%>%group_by(Year)%>%summarise(nm=n_distinct(Month))%>%ungroup()
as.data.frame(table(tab$nm))

tab2<-xweather%>%group_by(Year)%>%summarise(nm=unique(Month))%>%ungroup()
as.data.frame(table(tab2$nm))
# May-June-July mostly, but Apr to Aug present

climvar<-"tas"
cdat<-read.csv(here(paste("DATA/CHELSA_v2/monthly/",climvar,"_monthlyvalues_extracted_for_uRID_WP3.csv",sep="")))
cdat<-cdat%>%dplyr::select(-Longitude, -Latitude, -Stratum)

cnm<-colnames(cdat)
id4<-grep('tas_04_', cnm)
id5<-grep('tas_05_', cnm)
id6<-grep('tas_06_', cnm)
id7<-grep('tas_07_', cnm)
id8<-grep('tas_08_', cnm)
cdat<-cdat[,c(1,id4,id5,id6,id7,id8)]# extract months 
longcdat<- gather(cdat,key="key",value="value", -uRID)
longcdat<-longcdat%>%separate(key,c("var","month","year"))
longcdat$year<-as.integer(longcdat$year)

longcdat<-longcdat%>%group_by(uRID,year)%>%summarise(meanval=mean(value,na.rm=T))%>%ungroup()
#longcdat<-longcdat%>%group_by(uRID,year)%>%summarise(meanval=max(value,na.rm=T))%>%ungroup()


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

saveRDS(cdat,here("RESULTS/year_by_site_tas_avgAprtoAug.RDS"))

#=======================
prepare_tas_with_AprtoAug_Nyrs_threshold<-function(cdat,yr_threshold=40, climvar="tas"){
  
  if(yr_threshold==40){
    df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))
  }else{
    df<-read.csv(here(paste("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites_min",yr_threshold,"yr.csv",sep="")))
  }
  
  # These species are chosen
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
    
    # in celcius scale
    #mat_c<- -273.15 + (0.1*mat)
    #dmat_c<-pracma::detrend(mat_c)
    # We used the detrended data in further analysis, and the rank of it
    # so, I checked: copula::pobs(dmat[,1]) and copula::pobs(dmat_c[,1]) are exactly same
    # This shows the linear scale from kelvin to celcius transformation does not matter for ranked data
    
    d_allsite<-vector(mode="list",length=ncol(mat))
    names(d_allsite)<-colnames(mat)
    d_allsite_detrend<-d_allsite
    
    for(i in 1:length(d_allsite)){
      d_allsite[[i]]<-data.frame(Year=rownames(mat),Dat=mat[,i])
      d_allsite_detrend[[i]]<-data.frame(Year=rownames(dmat),Dat=dmat[,i])
    }
    
    saveRDS(mat,paste(resloc,"/year_by_site_",climvar,"_avgAprtoAug.RDS",sep=""))
    saveRDS(dmat,paste(resloc,"/year_by_site_",climvar,"_avgAprtoAug_detrended.RDS",sep=""))
    
    saveRDS(d_allsite,paste(resloc,"/",climvar,"_data_avgAprtoAug_selectedsitelist.RDS",sep=""))
    saveRDS(d_allsite_detrend,paste(resloc,"/",climvar,"_detrended_data_avgAprtoAug_selectedsitelist.RDS",sep=""))
    print(j)
  }
  
}
#===============
cdat<-readRDS(here("RESULTS/year_by_site_tas_avgAprtoAug.RDS"))
prepare_tas_with_AprtoAug_Nyrs_threshold(cdat=cdat, yr_threshold = 32,climvar="tas")
prepare_tas_with_AprtoAug_Nyrs_threshold(cdat=cdat, yr_threshold = 36,climvar="tas")
prepare_tas_with_AprtoAug_Nyrs_threshold(cdat=cdat, yr_threshold = 40,climvar="tas")

