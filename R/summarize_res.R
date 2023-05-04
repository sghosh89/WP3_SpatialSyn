library(tidyverse)
library(dplyr)
library(ggpubr)
library(here)

# this function keeps species set which have all finite non-zero value for flu_ab
# and finite value for spatial synchrony in climvar (e.g. pr, tas, tasmax, tasmin)

summarize_res<-function(chosen_rad){
  #============== read data ========================
  # the csv files you need for further analysis are:
  
  # spatial synchrony summary results (0-250 Km) for abundance, rainfall, temperature, maxT, minT
  
  
  df_ab<-read.csv(here(paste("RESULTS/summary_spat_syn_for_abund_",chosen_rad[1],
                             "_",chosen_rad[2],"km.csv",sep="")))
  df_pr<-read.csv(here(paste("RESULTS/summary_spat_syn_for_pr_",chosen_rad[1],
                             "_",chosen_rad[2],"km.csv",sep="")))
  df_tas<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tas_",chosen_rad[1],
                              "_",chosen_rad[2],"km.csv",sep="")))
  df_tasmax<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tasmax_",chosen_rad[1],
                                 "_",chosen_rad[2],"km.csv",sep="")))
  df_tasmin<-read.csv(here(paste("RESULTS/summary_spat_syn_for_tasmin_",chosen_rad[1],
                                 "_",chosen_rad[2],"km.csv",sep="")))
  
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
  
  # read tasmin summary
  df_tasmin<-df_tasmin%>%dplyr::select(AOU,tasmin_L=L,tasmin_U=U)%>%mutate(fLU_tasmin=(tasmin_L+tasmin_U)/(abs(tasmin_U)+tasmin_L))
  
  df<-cbind(df_ab$AOU,df_ab$fLU_ab,df_pr$fLU_pr,df_tas$fLU_tas,df_tasmax$fLU_tasmax,df_tasmin$fLU_tasmin)
  colnames(df)<-c("AOU","fLU_ab","fLU_pr","fLU_tas","fLU_tasmax","fLU_tasmin")
  df<-as.data.frame(df)
  
  
  # metadata: AOU code for given species, diet category and IUCN status
  df_spmeta<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
  df_spmeta<-df_spmeta%>%dplyr::select(AOU,Diet.5Cat,IUCN_status)
  
  df<-left_join(df,df_spmeta,by="AOU")# This is the dataframe we need to visualize
  
  
  id<-which(is.na(df$fLU_ab))
  df<-df[-id,]
  
  df<-na.omit(df) # still Na happens when no tail dep in any of the climate variables
  # e.g., 0-250km, fLU_pr shows NaN for AOU=7470
  
  df$tail<-ifelse(df$fLU_ab>0,"LT","UT")
  df$tail<-as.factor(df$tail)
  df$Diet.5Cat<-as.factor(df$Diet.5Cat)
  write.csv(df,here(paste("RESULTS/df_abund_climate_spatsyn_",
                          chosen_rad[1],"_",chosen_rad[2],"km.csv",sep="")),row.names = F)
  
  # you need this df dataframe file to plot abund_syn vs climate_syn
  
  #==================== plot across species ==============
  
  df<-read.csv(here(paste("RESULTS/df_abund_climate_spatsyn_",chosen_rad[1],
                          "_",chosen_rad[2],"km.csv",sep="")))
  ftmethod<-"lm"
  
  g1<-ggplot(data=df,aes(x=fLU_pr,y=fLU_ab,col=as.factor(tail)))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod)+xlab("Spatial synchrony, rainfall")+ylab("Spatial synchrony, abundance")+
    #facet_wrap(~Diet.5Cat)+
    geom_smooth(method=ftmethod, col="black")+
    theme_bw()+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "black", geom = "label")+ theme(legend.position="none")
  
  g2<-ggplot(data=df,aes(x=fLU_tas,y=fLU_ab,col=as.factor(tail)))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod)+xlab("Spatial synchrony, Temperature")+ylab("Spatial synchrony, abundance")+
    #facet_wrap(~Diet.5Cat)+
    geom_smooth(method=ftmethod, col="black")+
    theme_bw()+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "black", geom = "label")+ theme(legend.position="none")
  
  g3<-ggplot(data=df,aes(x=fLU_tasmax,y=fLU_ab,col=as.factor(tail)))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod)+xlab("Spatial synchrony, maximum T")+ylab("Spatial synchrony, abundance")+
    #facet_wrap(~Diet.5Cat)+
    geom_smooth(method=ftmethod, col="black")+
    theme_bw()+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "black", geom = "label")+ theme(legend.position="none")
  
  g4<-ggplot(data=df,aes(x=fLU_tasmin,y=fLU_ab,col=as.factor(tail)))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod)+xlab("Spatial synchrony, minimum T")+ylab("Spatial synchrony, abundance")+
    #facet_wrap(~Diet.5Cat)+
    geom_smooth(method=ftmethod, col="black")+
    theme_bw()+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "black", geom = "label")+ theme(legend.position="none")
  
  grid.arrange(g1, g2, g3, g4, ncol=2)
  
  #=============== plot groupwise: rare on top row, common on bottom row ===========
  ftmethod<-"lm"
  
  # ======== NOW WITH SPECIES WHICH ARE SYNCHRONOUSLY RARE (n=138) ===================
  # LT synchrony for abundance caused by high pr (+ve cor) and high temp across sites (-ve cor)
  
  df1<-df%>%filter(tail=="LT")
  g1<-ggplot(data=df1,aes(x=fLU_pr,y=fLU_ab))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod, aes(col="red"))+xlab("SpatSyn, rainfall")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+theme(legend.position = "none")+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "red", geom = "label")
  
  g2<-ggplot(data=df1,aes(x=fLU_tas,y=fLU_ab))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod,aes(col="red"))+xlab("SpatSyn, T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+theme(legend.position = "none")+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "red", geom = "label")
  
  g3<-ggplot(data=df1,aes(x=fLU_tasmax,y=fLU_ab))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod,aes(col="red"))+xlab("SpatSyn, max T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+theme(legend.position = "none")+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "red", geom = "label")
  
  g4<-ggplot(data=df1,aes(x=fLU_tasmin,y=fLU_ab))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod,aes(col="red"))+xlab("SpatSyn, min T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+theme(legend.position = "none")+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "red", geom = "label")
  
  # === NOW WITH SPECIES WHICH ARE SYNCHRONOUSLY COMMON (n=136, excluded AOU=7470; indep pr spatsyn for 2 sites) ===================
  # UT synchrony for abundance caused by high temp across sites (+ve cor)
  df2<-df%>%filter(tail=="UT")
  g5<-ggplot(data=df2,aes(x=fLU_pr,y=fLU_ab))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod)+xlab("SpatSyn, rainfall")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "steelblue3", geom = "label")
  
  g6<-ggplot(data=df2,aes(x=fLU_tas,y=fLU_ab))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod)+xlab("SpatSyn, T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "steelblue3", geom = "label")
  
  g7<-ggplot(data=df2,aes(x=fLU_tasmax,y=fLU_ab))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod)+xlab("SpatSyn, max T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "steelblue3", geom = "label")
  
  g8<-ggplot(data=df2,aes(x=fLU_tasmin,y=fLU_ab))+geom_point(alpha=0.1)+
    geom_smooth(method=ftmethod)+xlab("SpatSyn, min T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+stat_cor(aes(label = paste(after_stat(p.label), sep = "~")), color = "steelblue3", geom = "label")
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), color = "steelblue1", geom = "label")
  
  grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, nrow=2)
  
  #================ plot by diet group ================
  
  ftmethod<-"lm"
  
  # ======== NOW WITH SPECIES WHICH ARE SYNCHRONOUSLY RARE (n=138) ===================
  # LT synchrony for abundance caused by high pr (+ve cor) and high temp across sites (-ve cor)
  
  df1<-df%>%filter(tail=="LT")
  g1<-ggplot(data=df1,aes(x=fLU_pr,y=fLU_ab,col=Diet.5Cat))+geom_point(alpha=0.5)+
    geom_smooth(method=ftmethod,se=F)+xlab("SpatSyn, rainfall")+ylab("SpatSyn, abundance")+
    #facet_grid(~Diet.5Cat)+
    theme_bw()+ theme(legend.position="none")
  
  g2<-ggplot(data=df1,aes(x=fLU_tas,y=fLU_ab,col=Diet.5Cat))+geom_point(alpha=0.5)+
    geom_smooth(method=ftmethod,se=F)+xlab("SpatSyn, T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ theme(legend.position="none")
  
  g3<-ggplot(data=df1,aes(x=fLU_tasmax,y=fLU_ab,col=Diet.5Cat))+geom_point(alpha=0.5)+
    geom_smooth(method=ftmethod,se=F)+xlab("SpatSyn, max T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ theme(legend.position="none")
  
  g4<-ggplot(data=df1,aes(x=fLU_tasmin,y=fLU_ab,col=Diet.5Cat))+geom_point(alpha=0.5)+
    geom_smooth(method=ftmethod,se=F)+xlab("SpatSyn, min T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ theme(legend.position="none")
  
  # === NOW WITH SPECIES WHICH ARE SYNCHRONOUSLY COMMON (n=136, excluded AOU=7470; indep pr spatsyn for 2 sites) ===================
  # UT synchrony for abundance caused by high temp across sites (+ve cor)
  df2<-df%>%filter(tail=="UT")
  g5<-ggplot(data=df2,aes(x=fLU_pr,y=fLU_ab,col=Diet.5Cat))+geom_point(alpha=0.5)+
    geom_smooth(method=ftmethod,se=F)+xlab("SpatSyn, rainfall")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ theme(legend.position="none")
  
  g6<-ggplot(data=df2,aes(x=fLU_tas,y=fLU_ab,col=Diet.5Cat))+geom_point(alpha=0.5)+
    geom_smooth(method=ftmethod,se=F)+xlab("SpatSyn, T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ theme(legend.position="none")
  
  g7<-ggplot(data=df2,aes(x=fLU_tasmax,y=fLU_ab,col=Diet.5Cat))+geom_point(alpha=0.5)+
    geom_smooth(method=ftmethod,se=F)+xlab("SpatSyn, max T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ theme(legend.position="none")
  
  g8<-ggplot(data=df2,aes(x=fLU_tasmin,y=fLU_ab,col=Diet.5Cat))+geom_point(alpha=0.5)+
    geom_smooth(method=ftmethod,se=F)+xlab("SpatSyn, min T")+ylab("SpatSyn, abundance")+
    #facet_wrap(~Diet.5Cat)+
    theme_bw()+ theme(legend.position="none")
  
  g9<-ggplot(data=df2,aes(x=NA,y=NA,col=Diet.5Cat))+geom_point(alpha=0)+
    geom_smooth(method=ftmethod,se=F)+xlab("")+ylab("")+
    #facet_wrap(~Diet.5Cat)+
    theme_void()+ theme(legend.title=element_blank(),legend.position="top")
  
  grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, 
               layout_matrix = rbind(c(1, 2, 3, 4),
                                     c(5, 6, 7, 8),c(9, 9, 9,9)), nrow=3)
}

# input chosen radius 
#chosen_rad<-c(0,250)
#summarize_res(chosen_rad=chosen_rad)