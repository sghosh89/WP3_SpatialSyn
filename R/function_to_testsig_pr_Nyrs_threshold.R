#=========================================================================================
# We are testing the significance of corlmcoru for pr copula against normal copula
#==========================================================================================
rm(list=ls())

library(tidyverse)
library(dplyr)
library(here)
library(gridExtra)
library(ggpubr)
source(here("R/CorlCoru.R"))
source(here("R/copsurrog2d.R"))

testsig_pr_Nyrs_threshold<-function(numsurr=1000, nbin=4, chosen_rad=c(0,250), yr_threshold=40){
  
  if(yr_threshold==40){
    df_pr<-read.csv(here(paste("RESULTS/summary_spat_syn_for_pr_",chosen_rad[1],
                               "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
    df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4.csv"))
  }else{
    df_pr<-read.csv(here(paste("RESULTS/summary_spat_syn_for_pr_",chosen_rad[1],
                                "_",chosen_rad[2],"km_nbin_",nbin,"_min",yr_threshold,"yr.csv",sep="")))
    df<-read.csv(here(paste("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_",nbin,"_min",yr_threshold,"yr.csv",sep="")))
  }
  
  
  df_pr<-df_pr%>%dplyr::select(AOU,nL,nU,L,U)
  
  df<-df%>%dplyr::select(AOU,fLU_pr)
  df<-left_join(df,df_pr,by="AOU") # this is summarized value within 0-250Km distance
  
  #=============
  for(spi in 1:nrow(df)){
    
    set.seed(seed=101)
    
    #spi<-2
    
    givenAOU<-df$AOU[spi]
    
    if(yr_threshold==40){
      inputdistloc<-here(paste("RESULTS/AOU_", givenAOU,sep=""))
      inputresloc<-here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn_nbin_",nbin,sep=""))
    }else{
      inputdistloc<-here(paste("RESULTS/yr_threshold_",yr_threshold,"/AOU_", givenAOU,sep=""))
      inputresloc<-here(paste("RESULTS/yr_threshold_",yr_threshold,"/AOU_", givenAOU,"/abundance_spatsyn_nbin_",nbin,sep=""))
    }
    
    sigabund<-readRDS(here(paste(inputresloc,"/corlmcoru_sigres.RDS",sep="")))
    sigabund$row_col<-paste(sigabund$row,sigabund$col,sep="_")
    
    #------------- distance category to choose -------------------
    distm<-readRDS(here(paste(inputdistloc,"/distm_sel.RDS",sep="")))
    # now only select sites which are within the interaction radius
    badid_distm_sel<-which(distm<chosen_rad[1] | distm>chosen_rad[2],arr.ind = T)
    distm_good<-distm
    distm_good[badid_distm_sel]<-NA
    
    #--------------corl - coru within that distance category showing +ve correlation----------
    nps<-readRDS(here(paste(inputdistloc,"/pr_spatsyn_nbin_",nbin,
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
    
    dat_detrend<-readRDS(here(paste(inputdistloc,"/year_by_site_pr_detrended.RDS",sep="")))
    
    #i<-1
    for(i in 1:nrow(dfid)){
      dr<-dfid$row[i]
      dc<-dfid$col[i]
      rawdat<-dat_detrend[,c(dr,dc)]
      rawdat<-na.omit(rawdat) # omitting NA value in 2019 for pr
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
    
    saveRDS(dfid,here(paste(inputdistloc,"/pr_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    print(spi)
  }
}
#=============================
# Now call the function
testsig_pr_Nyrs_threshold(numsurr = 1000, nbin=4, chosen_rad = c(0,250), yr_threshold = 32)
testsig_pr_Nyrs_threshold(numsurr = 1000, nbin=4, chosen_rad = c(0,250), yr_threshold = 36)
testsig_pr_Nyrs_threshold(numsurr = 1000, nbin=4, chosen_rad = c(0,250), yr_threshold = 40)













