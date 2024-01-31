# Make a summary table for each species 
# with sig (75%CI) taildep in abundance and climate (tasmax) for each site pair distance

rm(list=ls())
library(here)
library(tidyverse)
nbin<-4
chosen_rad<-c(0,250) # within this distance category

df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4.csv"))
df<-df%>%dplyr::select(AOU,fLU_ab,fLU_tasmax,fLU_tasmax_avgMaytoJuly) # before considering significance

# Now consider tail-dep significant result

distance_sigtaildep_abund_tasmax<-c()
for(i in 1:nrow(df)){
  
  givenAOU<-df$AOU[i]
  sigabund<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
  sigabund$row_col<-paste(sigabund$row,sigabund$col,sep="_")
  sigabund<-sigabund%>%dplyr::select(row_col,corlmcoru_ab=corlmcoru_actual,
                                     sig75ab=sig75,sig95ab=sig95)
  idgs<-readRDS(paste("RESULTS/AOU_",givenAOU,"/abundance_spatsyn_nbin_",nbin,"/abund_table_for_sitepair_within_",chosen_rad[1],"_",chosen_rad[2],"_km_nbin_",nbin,".RDS",sep=""))
  
  idgs<-idgs%>%dplyr::mutate(row_col=paste(row,col,sep="_"))
  idgs<-idgs%>%dplyr::select(row_col,dist.KM)
  sigabund<-left_join(sigabund,idgs,by="row_col")
  
  #------ for climate: tasmax avg across all 12 months --------
  sigclim<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/tasmax_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
  sigclim<-sigclim%>%dplyr::select(row_col,corlmcoru_tasmax=corlmcoru_actual,
                                     sig75tasmax=sig75,sig95tasmax=sig95)
  
  sigabund<-left_join(sigabund,sigclim,by="row_col")
  
  #------ for climate: tasmax avg across 3 months: May to July --------
  sigclim3<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/tasmax_avgMaytoJuly_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
  sigclim3<-sigclim3%>%dplyr::select(row_col,corlmcoru_tasmax3=corlmcoru_actual,
                                   sig75tasmax3=sig75,sig95tasmax3=sig95)
  sigabund<-left_join(sigabund,sigclim3,by="row_col")
  
  sigabund$AOU<-givenAOU
  
  sigabund<-sigabund%>%filter(sig75ab==1)
  distance_sigtaildep_abund_tasmax<-rbind(distance_sigtaildep_abund_tasmax,sigabund)
}
distance_sigtaildep_abund_tasmax<-distance_sigtaildep_abund_tasmax%>%
                                  dplyr::select(AOU,row_col,dist.KM,
                                                corlmcoru_ab,
                                                corlmcoru_tasmax,
                                                corlmcoru_tasmax3,
                                                sig75ab,
                                                sig75tasmax,
                                                sig75tasmax3)

write.csv(distance_sigtaildep_abund_tasmax,here("RESULTS/distance_sigtaildep_abund_tasmax.csv"), row.names = F)
#==================================================
rm(list=ls())
library(here)
library(tidyverse)
dat<-read.csv(here("RESULTS/distance_sigtaildep_abund_tasmax.csv"))

# split into two category: <=100 km and 100-250 km
df100<-dat%>%filter(dist.KM<100)
df<-data.frame(AOU=unique(df100$AOU))

# values: total sig within 0-100 Km distance
df$fab.sig<-NA
df$ftasmax.sig<-NA

for(i in 1:nrow(df)){
  
  givenAOU<-df$AOU[i]
  tempo<-df100%>%filter(AOU%in%givenAOU)
  df$fab.sig[i]<-sum(tempo$corlmcoru_ab)/sum(abs(tempo$corlmcoru_ab))
  
  tempoclim<-tempo%>%filter(sig75tasmax==1)
  df$ftasmax.sig[i]<-sum(tempoclim$corlmcoru_ab)/sum(abs(tempoclim$corlmcoru_ab))
}

df$tail<-ifelse(df$fab.sig>0,"LT","UT")

# <100 Km. category
table(df$tail)# 44 (21 LT, 23 UT) sp. out of 59 species below 100 Km distance

# 100-250 Km. category
#table(df$tail)# 50 (29 LT, 21 UT)sp. out of 59 species fall within 100-250 Km distance

# <100 Km. category
df<-na.omit(df) # 24 sp. left after we delete NaN values of climate
table(df$tail) # 12 LT, 12 UT

# 100-250 Km. category
#df<-na.omit(df) # 30 sp. left after we delete NaN values of climate
#table(df$tail) # 20 LT, 10 UT








