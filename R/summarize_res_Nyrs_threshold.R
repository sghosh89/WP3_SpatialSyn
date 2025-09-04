library(tidyverse)
library(dplyr)
library(ggpubr)
library(here)
library(gridExtra)

# this function keeps species set which have all finite non-zero value for flu_ab
# and finite value for spatial synchrony in climvar (e.g. pr, tasmax)

summarize_res_Nyrs_threshold<-function(yr_threshold=40,chosen_rad,nbin){
  #============== read data ========================
  # the csv files you need for further analysis are:
  
  # spatial synchrony summary results (0-250 Km) for abundance, Precipitation, temperature
  
  if(yr_threshold==40){
    df_ab<-read.csv(here(paste("RESULTS/summary_spat_syn_for_abund_",chosen_rad[1],
                               "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
    
    df_pr<-read.csv(here(paste("RESULTS/summary_spat_syn_for_pr_",chosen_rad[1],
                               "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
    df_prAprtoAugavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_pr_avgAprtoAug_",chosen_rad[1],
                                          "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
    
    
    df_tas<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_",chosen_rad[1],
                                "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
    df_tasAprtoAugavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_avgAprtoAug_",chosen_rad[1],
                                           "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
    #df_tasAprtoJulyavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_avgAprtoJuly_",chosen_rad[1],
    #                                       "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
    #df_tasMaytoJulyavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_avgMaytoJuly_",chosen_rad[1],
    #                                        "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
    
  }else{
    
    df_ab<-read.csv(here(paste("RESULTS/summary_spat_syn_for_abund_",chosen_rad[1],
                               "_",chosen_rad[2],"km_nbin_",nbin,"_min",yr_threshold,"yr.csv",sep="")))
    
    df_pr<-read.csv(here(paste("RESULTS/summary_spat_syn_for_pr_",chosen_rad[1],
                               "_",chosen_rad[2],"km_nbin_",nbin,"_min",yr_threshold,"yr.csv",sep="")))
    df_prAprtoAugavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_pr_avgAprtoAug_",chosen_rad[1],
                                          "_",chosen_rad[2],"km_nbin_",nbin,"_min",yr_threshold,"yr.csv",sep="")))
    
    
    df_tas<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_",chosen_rad[1],
                                "_",chosen_rad[2],"km_nbin_",nbin,"_min",yr_threshold,"yr.csv",sep="")))
    df_tasAprtoAugavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_avgAprtoAug_",chosen_rad[1],
                                           "_",chosen_rad[2],"km_nbin_",nbin,"_min",yr_threshold,"yr.csv",sep="")))
  }
  
  #======================== make combo data ===========================
  # read abundance summary
  df_ab<-df_ab%>%dplyr::select(AOU,ab_L=L,ab_U=U)%>%mutate(fLU_ab=(ab_L+ab_U)/(abs(ab_U)+ab_L))
  
  # read pr summary
  df_pr<-df_pr%>%dplyr::select(AOU,pr_L=L,pr_U=U)%>%mutate(fLU_pr=(pr_L+pr_U)/(abs(pr_U)+pr_L))
  
  # read prs 5 months avg summary
  df_prAprtoAugavg<-df_prAprtoAugavg%>%dplyr::select(AOU,pr5_L=L,pr5_U=U)%>%mutate(fLU_pr5=(pr5_L+pr5_U)/(abs(pr5_U)+pr5_L))
  
  # read tas summary
  df_tas<-df_tas%>%dplyr::select(AOU,tas_L=L,tas_U=U)%>%mutate(fLU_tas=(tas_L+tas_U)/(abs(tas_U)+tas_L))
  
  # read tas 5 months avg summary
  df_tasAprtoAugavg<-df_tasAprtoAugavg%>%dplyr::select(AOU,tas5_L=L,tas5_U=U)%>%mutate(fLU_tas5=(tas5_L+tas5_U)/(abs(tas5_U)+tas5_L))
  
  df<-cbind(df_ab$AOU,df_ab$fLU_ab,
            df_pr$fLU_pr,
            df_tas$fLU_tas,
            df_prAprtoAugavg$fLU_pr5,
            df_tasAprtoAugavg$fLU_tas5)
  colnames(df)<-c("AOU","fLU_ab","fLU_pr","fLU_tas",
                  "fLU_pr_avgAprtoAug",
                  "fLU_tas_avgAprtoAug")
  df<-as.data.frame(df)
  
  
  # metadata: AOU code for given species, diet category and IUCN status
  #df_spmeta<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
  #df_spmeta<-df_spmeta%>%dplyr::select(AOU,Diet.5Cat,IUCN_status)
  
  #df<-left_join(df,df_spmeta,by="AOU")# This is the dataframe we need to visualize
  
  
  id<-which(is.na(df$fLU_ab))
  df<-df[-id,]
  
  df<-na.omit(df) # still NA happens when no tail dep in any of the climate variables
  # e.g., 0-250km, fLU_pr shows NaN for AOU=7470 for 40yrs. 
  
  df$tail<-ifelse(df$fLU_ab>0,"LT","UT")
  df$tail<-as.factor(df$tail)
  #df$Diet.5Cat<-as.factor(df$Diet.5Cat)
  
  if(yr_threshold==40){
    write.csv(df,here(paste("RESULTS/df_abund_climate_spatsyn_",
                            chosen_rad[1],"_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")),row.names = F)
  }else{
    write.csv(df,here(paste("RESULTS/df_abund_climate_spatsyn_",
                            chosen_rad[1],"_",chosen_rad[2],"km_nbin_",nbin,"_min",yr_threshold,"yr.csv",sep="")),row.names = F)
  }
 
  # you need this df dataframe file to plot abund_syn vs climate_syn
  return(df)
}

# input chosen radius 
chosen_rad<-c(0,250)
df32<-summarize_res_Nyrs_threshold(yr_threshold=32,chosen_rad=chosen_rad,nbin=4)#236sp
df36<-summarize_res_Nyrs_threshold(yr_threshold=36,chosen_rad=chosen_rad,nbin=4)#168sp
df40<-summarize_res_Nyrs_threshold(yr_threshold=40,chosen_rad=chosen_rad,nbin=4)#78sp


chosen_rad<-c(0,100)
df40a<-summarize_res_Nyrs_threshold(yr_threshold=40,chosen_rad=chosen_rad,nbin=4)#62sp

chosen_rad<-c(100,250)
df40b<-summarize_res_Nyrs_threshold(yr_threshold=40,chosen_rad=chosen_rad,nbin=4)#71sp

#df36$AOU%in%df32$AOU # this means all 168 species in df36 data are a subset of 236 species in df32 data
#df40$AOU%in%df32$AOU # this means all 78 species in df36 data are a subset of 236 species in df32 data
# therefore we need trait info for those 236 species


