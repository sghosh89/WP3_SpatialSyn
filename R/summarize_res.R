library(tidyverse)
library(dplyr)
library(ggpubr)
library(here)
library(gridExtra)

# this function keeps species set which have all finite non-zero value for flu_ab
# and finite value for spatial synchrony in climvar (e.g. pr, tasmax)

summarize_res<-function(chosen_rad,nbin){
  #============== read data ========================
  # the csv files you need for further analysis are:
  
  # spatial synchrony summary results (0-250 Km) for abundance, Precipitation, temperature
  
  
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
  
  # read tas 4 months avg summary
  #df_tasAprtoJulyavg<-df_tasAprtoJulyavg%>%dplyr::select(AOU,tas4_L=L,tas4_U=U)%>%mutate(fLU_tas4=(tas4_L+tas4_U)/(abs(tas4_U)+tas4_L))
  
  # read tas 3 months avg summary
  #df_tasMaytoJulyavg<-df_tasMaytoJulyavg%>%dplyr::select(AOU,tas3_L=L,tas3_U=U)%>%mutate(fLU_tas3=(tas3_L+tas3_U)/(abs(tas3_U)+tas3_L))
  
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
  df_spmeta<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
  df_spmeta<-df_spmeta%>%dplyr::select(AOU,Diet.5Cat,IUCN_status)
  
  df<-left_join(df,df_spmeta,by="AOU")# This is the dataframe we need to visualize
  
  
  id<-which(is.na(df$fLU_ab))
  df<-df[-id,]
  
  df<-na.omit(df) # still NA happens when no tail dep in any of the climate variables
  # e.g., 0-250km, fLU_pr shows NaN for AOU=7470
  
  df$tail<-ifelse(df$fLU_ab>0,"LT","UT")
  df$tail<-as.factor(df$tail)
  df$Diet.5Cat<-as.factor(df$Diet.5Cat)
  write.csv(df,here(paste("RESULTS/df_abund_climate_spatsyn_",
                          chosen_rad[1],"_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")),row.names = F)
  
  # you need this df dataframe file to plot abund_syn vs climate_syn
  return(df)
}

# input chosen radius 
chosen_rad<-c(0,250)
df4<-summarize_res(chosen_rad=chosen_rad,nbin=4)




