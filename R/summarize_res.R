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
  df_tas<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_",chosen_rad[1],
                              "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
  
  df_tasMaytoJulyavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_avgMaytoJuly_",chosen_rad[1],
                                "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
  df_tasAprtoAugavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_avgAprtoAug_",chosen_rad[1],
                                            "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
  df_tasAprtoJulyavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_avgAprtoJuly_",chosen_rad[1],
                                         "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
  
  
  df_tasmax<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tasmax_",chosen_rad[1],
                                 "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
  #df_tasmaxMaytoJulyavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tasmax_avgMaytoJuly_",chosen_rad[1],
  #                              "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
  #df_tasmaxAprtoAugavg<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tasmax_avgAprtoAug_",chosen_rad[1],
  #                                           "_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
                                             
  #======================== make combo data ===========================
  # read abundance summary
  #df_ab<-read.csv(here(paste("RESULTS/summary_spat_syn_for_abund_",chosen_rad[1],
  #                           "_",chosen_rad[2],"km.csv",sep="")))
  df_ab<-df_ab%>%dplyr::select(AOU,ab_L=L,ab_U=U)%>%mutate(fLU_ab=(ab_L+ab_U)/(abs(ab_U)+ab_L))
  
  # read pr summary
  df_pr<-df_pr%>%dplyr::select(AOU,pr_L=L,pr_U=U)%>%mutate(fLU_pr=(pr_L+pr_U)/(abs(pr_U)+pr_L))
  
  # read tas summary
  df_tas<-df_tas%>%dplyr::select(AOU,tas_L=L,tas_U=U)%>%mutate(fLU_tas=(tas_L+tas_U)/(abs(tas_U)+tas_L))
  
  # read tasmax summary
  df_tasmax<-df_tasmax%>%dplyr::select(AOU,tasmax_L=L,tasmax_U=U)%>%mutate(fLU_tasmax=(tasmax_L+tasmax_U)/(abs(tasmax_U)+tasmax_L))
  
  # read tas 5 months avg summary
  df_tasAprtoAugavg<-df_tasAprtoAugavg%>%dplyr::select(AOU,tas5_L=L,tas5_U=U)%>%mutate(fLU_tas5=(tas5_L+tas5_U)/(abs(tas5_U)+tas5_L))
  
  # read tas 4 months avg summary
  df_tasAprtoJulyavg<-df_tasAprtoJulyavg%>%dplyr::select(AOU,tas4_L=L,tas4_U=U)%>%mutate(fLU_tas4=(tas4_L+tas4_U)/(abs(tas4_U)+tas4_L))
  
  # read tas 3 months avg summary
  df_tasMaytoJulyavg<-df_tasMaytoJulyavg%>%dplyr::select(AOU,tas3_L=L,tas3_U=U)%>%mutate(fLU_tas3=(tas3_L+tas3_U)/(abs(tas3_U)+tas3_L))
  
  df<-cbind(df_ab$AOU,df_ab$fLU_ab,df_pr$fLU_pr,df_tas$fLU_tas,
            df_tasmax$fLU_tasmax, 
            df_tasAprtoAugavg$fLU_tas5,
            df_tasAprtoJulyavg$fLU_tas4,
            df_tasMaytoJulyavg$fLU_tas3)
  colnames(df)<-c("AOU","fLU_ab","fLU_pr","fLU_tas","fLU_tasmax",
                  "fLU_tas_avgAprtoAug",
                  "fLU_tas_avgAprtoJuly",
                  "fLU_tas_avgMaytoJuly")
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




# run later below code
#==================== plot across species ==============

nbin<-4
df<-df4

df$tail<-factor(df$tail)
ftmethod<-"lm"
  
g1<-ggplot(data=df,aes(x=fLU_pr,y=fLU_ab), add = "reg.line")+
  geom_point(pch=21, col="white")+
    #geom_smooth(method=ftmethod,linetype="dashed")+
  xlab("Tail-dep. spatial synchrony, Precipitation")+
  ylab("Tail-dep. spatial synchrony, Abundance")+
    #facet_wrap(~Diet.5Cat)+
    geom_smooth(method=ftmethod, col="black")+
    theme_bw()+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = -0.7, label.y = 0.8,col="black")+
    stat_regline_equation(label.x = -0.7, label.y = 0.9)+
    geom_point(data=df, aes(x=fLU_pr,y=fLU_ab, col=tail), alpha=0.3)+
    theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g2<-ggplot(data=df,aes(x=fLU_tasmax,y=fLU_ab), add = "reg.line")+
  geom_point(pch=21, col="white")+
  #geom_smooth(method=ftmethod,linetype="dashed")+
  xlab("Tail-dep. spatial synchrony, Max. T")+
  ylab("Tail-dep. spatial synchrony, Abundance")+
  #facet_wrap(~Diet.5Cat)+
  geom_smooth(method=ftmethod, col="black")+
  theme_bw()+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = -0.7, label.y = 0.8,col="black")+
  stat_regline_equation(label.x = -0.7, label.y = 0.9)+
  geom_point(data=df, aes(x=fLU_tas,y=fLU_ab, col=tail), alpha=0.3)+
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  grid.arrange(g1, g2, ncol=2)
  # saving the plot as pdf
  pdf(here(paste("RESULTS/df_abund_climate_spatsyn_altogether_",
                 chosen_rad[1],"_",chosen_rad[2],"km_nbin_",nbin,".pdf",sep="")), width = 10, height = 4)
  grid.arrange(g1, g2, ncol=2)
  dev.off()
#=============== plot groupwise: rare on top row, common on bottom row ===========
  ftmethod<-"lm"
  
# ======== NOW WITH SPECIES WHICH ARE SYNCHRONOUSLY RARE (n=138) ===================
  # LT synchrony for abundance caused by high pr (+ve cor) and high temp across sites (-ve cor)
  
  df1<-df%>%filter(tail=="LT")
  g1<-ggplot(data=df1,aes(x=fLU_pr,y=fLU_ab), add = "reg.line"
             )+geom_point(alpha=0.3, col = "red")+
    geom_smooth(method=ftmethod, aes(col="red"))+
    xlab("Tail-dep. spatial synchrony, Precipitation")+
    ylab("Tail-dep. spatial synchrony, Abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ 
    theme(legend.position="none",panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
             label.x = -0.7, label.y = 0.8,col="black")+
    stat_regline_equation(label.x = -0.7, label.y = 0.9)
  
  #g1
  
  
  g2<-ggplot(data=df1,aes(x=fLU_tasmax,y=fLU_ab), add = "reg.line")+geom_point(alpha=0.3, col = "red")+
    geom_smooth(method=ftmethod,aes(col="red"))+
    xlab("Tail-dep. spatial synchrony, Max. T")+
    ylab("Tail-dep. spatial synchrony, Abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ 
    theme(legend.position="none",panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
             label.x = -0.7, label.y = 0.8,col="black")+
    stat_regline_equation(label.x = -0.7, label.y = 0.9)
  #g2
  
# === NOW WITH SPECIES WHICH ARE SYNCHRONOUSLY COMMON ===================
  # UT synchrony for abundance caused by high temp across sites (+ve cor)
  df2<-df%>%filter(tail=="UT")
  g5<-ggplot(data=df2,aes(x=fLU_pr,y=fLU_ab))+geom_point(alpha=0.3, col = "darkturquoise")+
    geom_smooth(method=ftmethod, col = "darkturquoise")+
    xlab("Tail-dep. spatial synchrony, Precipitation")+
    ylab("Tail-dep. spatial synchrony, Abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ 
    theme(legend.position="none",panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
             label.x = -0.7, label.y = -0.9,col="black")+
    stat_regline_equation(label.x = -0.7, label.y = -0.8)
  g5
  
  
  g6<-ggplot(data=df2,aes(x=fLU_tasmax,y=fLU_ab))+geom_point(alpha=0.3, col = "darkturquoise")+
    geom_smooth(method=ftmethod, col = "darkturquoise")+
    xlab("Tail-dep. spatial synchrony, Max. T")+
    ylab("Tail-dep. spatial synchrony, Abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ 
    theme(legend.position="none",panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
             label.x = -0.7, label.y = -0.9,col="black")+
    stat_regline_equation(label.x = -0.7, label.y = -0.8)
  g6
  
  grid.arrange(g1, g2, g5, g6, nrow=2)
  
  # saving the plot as pdf
  pdf(here(paste("RESULTS/df_abund_climate_spatsyn_groupwise_",
                 chosen_rad[1],"_",chosen_rad[2],"km_nbin_",nbin,".pdf",sep="")), 
      width = 8, height = 7)
  
  grid.arrange(g1, g2, g5, g6, nrow=2)
  dev.off()
