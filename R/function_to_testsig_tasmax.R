#=========================================================================================
# We are testing the significance of corlmcoru for tasmax copula against normal copula
#==========================================================================================
rm(list=ls())

library(tidyverse)
library(dplyr)
library(here)
library(gridExtra)
library(ggpubr)
source(here("R/CorlCoru.R"))
source(here("R/copsurrog2d.R"))

numsurr<-1000

nbin<-4
chosen_rad<-c(0,250) # within this distance category

df_tasmax<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tasmax_",chosen_rad[1],
                           "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
df_tasmax<-df_tasmax%>%dplyr::select(AOU,nL,nU,L,U)

df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4.csv"))
df<-df%>%dplyr::select(AOU,fLU_tasmax)
df<-left_join(df,df_tasmax,by="AOU") # this is summarized value within 0-250Km distance

#=============
for(spi in 1:nrow(df)){
  
  set.seed(seed=101)
  
  #spi<-2
  
  givenAOU<-df$AOU[spi]
  sigabund<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
  sigabund$row_col<-paste(sigabund$row,sigabund$col,sep="_")
  
  #------------- distance category to choose -------------------
  distm<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/distm_sel.RDS",sep="")))
  # now only select sites which are within the interaction radius
  badid_distm_sel<-which(distm<chosen_rad[1] | distm>chosen_rad[2],arr.ind = T)
  distm_good<-distm
  distm_good[badid_distm_sel]<-NA
  
  #--------------corl - coru within that distance category showing +ve correlation----------
  
  nps<-readRDS(here(paste("RESULTS/AOU_",givenAOU,"/tasmax_spatsyn_nbin_",nbin,
                          "/NonParamStat.RDS",sep="")))
  
  corlmcoru_good<- nps$Corl - nps$Coru
  corlmcoru_good[badid_distm_sel]<-NA # exclude outside of chosen radius
  
  posnImat_old<-nps$posnI
  posnImat_new<-posnImat_old
  posnImat_new[badid_distm_sel]<-NA
  posnI_new<-which(posnImat_new==1,arr.ind = T)
  
  posnNmat_old<-nps$posnN
  posnNmat_new<-posnNmat_old
  posnNmat_new[badid_distm_sel]<-NA
  posnN_new<-which(posnNmat_new==1,arr.ind = T)
  
  corlmcoru_good[posnI_new]<-NA
  corlmcoru_good[posnN_new]<-NA
  
  dfid<-which(!is.na(corlmcoru_good),arr.ind=T)# this could have multiple rows , so put a for loop
  dfid<-as.data.frame(dfid)
  dfid$corlmcoru_actual<-NA
  dfid$sig75<-NA
  dfid$sig95<-NA
  dfid$row_col<-paste(dfid$row,dfid$col,sep="_")
  idm<-which(dfid$row_col%in%sigabund$row_col) # only consider cells that matches with abundance spat syn
  dfid<-dfid[idm,]
  #-------------------------------------------
  
  dat_detrend<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/year_by_site_tasmax_detrended.RDS",sep="")))
  
  #i<-1
  for(i in 1:nrow(dfid)){
    dr<-dfid$row[i]
    dc<-dfid$col[i]
    rawdat<-dat_detrend[,c(dr,dc)]
    rdat<-copula::pobs(rawdat)
    
    corlmcoru_actual<-corlmcoru_good[dr,dc]
    
    # first generate Normal surrogs array
    tcop<-copula::normalCopula(0.5,2)
    
    sarray <- copsurrog2d(m = rdat, targetcop = tcop, corpres = "spearman", numsurr = numsurr)
    
    corl_s <- coru_s<- corlmcoru_s<- c()
    for(ir in 1:numsurr){
      sdat <- sarray[ , , ir]
      call_on_sdat <- CorlCoru(vi=sdat[,1],vj=sdat[,2],nbin=nbin)
      corl_s <- c(corl_s, call_on_sdat[1])
      coru_s <- c(coru_s, call_on_sdat[2])
      tempo<-call_on_sdat[1]-call_on_sdat[2]
      corlmcoru_s<-c(corlmcoru_s, tempo)
    }
    
    #hist(corlmcoru_s, col="white",1000)
    #abline(v=corlmcoru_actual,col="red")
    #abline(v=ci95[1],col="blue")
    #abline(v=ci95[2],col="blue")
    
    dfid$corlmcoru_actual[i]<-corlmcoru_actual
    
    ci75<-quantile(corlmcoru_s,probs=c(0.125,0.875))
    nsig75<-(corlmcoru_actual>as.numeric(ci75[1]) & corlmcoru_actual<as.numeric(ci75[2]))# value within ci means not sig
    dfid$sig75[i]<- 1-nsig75
    
    ci95<-quantile(corlmcoru_s,probs=c(0.025,0.975))
    nsig95<-(corlmcoru_actual>as.numeric(ci95[1]) & corlmcoru_actual<as.numeric(ci95[2]))# value within ci means not sig
    dfid$sig95[i]<- 1-nsig95
    
  }
  
  myind<-as.matrix(dfid[,c(1,2)])
  mydist<-distm_good[myind]
  dfid$dist.KM<-mydist
  
  saveRDS(dfid,here(paste("RESULTS/AOU_", givenAOU,"/tasmax_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
  print(spi)
}
#=============================


















#=================== Now, sum significant taildep. within 0-250 km distance =======
df$Lsig75<-NA
df$Usig75<-NA
df$Lsig95<-NA
df$Usig95<-NA
for(i in 1:nrow(df)){
  
  givenAOU<-df$AOU[i]
  tempo<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/tasmax_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
  
  tempo75<-tempo%>%filter(sig75==1)
  
  if(nrow(tempo75)==0){
    Lsig75<-Usig75<-0
  }else{
    Lsig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual>0)])
    Usig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual<0)])
  }
  
  tempo95<-tempo%>%filter(sig95==1)
  
  if(nrow(tempo95)==0){
    Lsig95<-Usig95<-0
  }else{
    Lsig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual>0)])
    Usig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual<0)])
  }
  
  df$Lsig75[i]<-Lsig75
  df$Usig75[i]<-Usig75
  df$Lsig95[i]<-Lsig95
  df$Usig95[i]<-Usig95
}

saveRDS(df,here(paste("RESULTS/tasmax_spatsyn_nbin_",nbin,"_corlmcoru_sigres_summary.RDS",sep="")))


