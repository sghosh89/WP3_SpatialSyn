# plot spatial syn (CorL and CorU) for species against distance category
#rm(list=ls())
library(here)
library(dplyr)
library(tidyverse)
library(gridExtra)

call_spat_syn_vs_distance_plot_Nyrs_threshold<-function(yr_threshold=40){
  
  if(yr_threshold==40){
    df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))
  }else{
    df<-read.csv(here(paste("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites_min",yr_threshold,"yr.csv",sep="")))
  }
  
  dflong<-c()
  
  for(i in 1:nrow(df)){
    
    givenAOU<-df$AOU[i]
    
    if(yr_threshold==40){
      distm<-readRDS(here(paste("RESULTS/AOU_",givenAOU,"/distm_sel.RDS",sep="")))
      tempo<-readRDS(here(paste("RESULTS/AOU_",givenAOU,"/abundance_spatsyn_nbin_4/NonParamStat.RDS",sep="")))
    }else{
      distm<-readRDS(here(paste("RESULTS/yr_threshold_",yr_threshold,"/AOU_",givenAOU,"/distm_sel.RDS",sep="")))
      tempo<-readRDS(here(paste("RESULTS/yr_threshold_",yr_threshold,"/AOU_",givenAOU,"/abundance_spatsyn_nbin_4/NonParamStat.RDS",sep="")))
    }
    
    
    indI<-which(tempo$posnI==1,arr.ind = T)
    indN<-which(tempo$posnN==1,arr.ind = T)
    
    cl<-tempo$Corl
    cl[indI]<-NA
    cl[indN]<-NA
    
    cu<-tempo$Coru
    cu[indI]<-NA
    cu[indN]<-NA
    
    spear<-tempo$corval #  spatial syn = +ve/-ve both
    spear[indI]<-NA
    #spear[indN]<-NA
    #plot(distm,cu)
    
    dtab<-data.frame(d=as.vector(distm),cl=as.vector(cl),cu=as.vector(cu), spear=as.vector(spear))
    dtab<-na.omit(dtab)
    
    if(nrow(dtab!=0)){
      dtab$AOU<-givenAOU
      dtab<-dtab%>%arrange(d)
      dflong<-rbind(dflong,dtab)
    }
    print(i)
  }
  
  # now plot
  range(dflong$d)
  
  dflong$dc<-cut(dflong$d, breaks=c(0,25,75,150,250,400,700,1000,1500,2500,6000))
  
  #as.data.frame(table(dflong$dc))
  
  dfextra<-dflong%>%group_by(dc)%>%summarise(n=n_distinct(AOU))%>%ungroup()
  #=======================
  
  dfs<-dflong%>%group_by(dc)%>%summarise(cl=mean(cl),cu=mean(cu),spear=mean(spear))%>%ungroup()
  
  dg1<-dfs%>%dplyr::select(dc,value=cl)%>%mutate(type="Corl")
  dg2<-dfs%>%dplyr::select(dc,value=cu)%>%mutate(type="Coru")
  dg3<-dfs%>%dplyr::select(dc,value=spear)%>%mutate(type="Spearman Cor")
  dg<-rbind(dg1,dg2,dg3)
  
  g1<-ggplot(dg, aes(x=dc, y=value, col=type)) +
    geom_point()+
    #geom_line()+
    xlab("Pairwise-distance category, Km")+ylab("Spatial synchrony")+
    theme_bw()+theme(legend.position = c(0.7,0.8))
  #=======================
  
  df2<-dflong%>%dplyr::select(AOU,dc,spear)
  
  # first average over all site-pair interactions for a given distance category and for a given species
  dff<-df2%>%group_by(dc,AOU)%>%summarise(mn.spear=mean(spear))%>%ungroup()
  
  
  library(PupillometryR)
  dff$AOU<-as.factor(dff$AOU)
  
  dff$numericdc<-as.numeric(dff$dc)
  dff2<-dff
  
  g2ex<-ggplot(dff, aes(x=dc, y=mn.spear)) + 
    #geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 1) + 
    geom_point(aes(y = mn.spear), 
               position = position_jitter(width = .15), size = 1, alpha = 0.3) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.2)+ 
    #geom_smooth(method="loess")+
    theme_bw()+
    theme(
      #panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
      panel.background=element_rect(fill="white", colour="white"), 
      legend.position="none",text=element_text(size=20)
    )+
    geom_smooth(data=dff2,aes(x=numericdc, y=mn.spear),method="loess")+
    ylab("Spatial synchrony in abundance\n averaged across site-pairs")+
    xlab("Between-sites pairwise distance category, Km") 
  
  
  
  return(list(tab=dfextra,gp=g2ex))
  
}

yr_threshold<-32
p32<-call_spat_syn_vs_distance_plot_Nyrs_threshold(yr_threshold = yr_threshold)

yr_threshold<-36
p36<-call_spat_syn_vs_distance_plot_Nyrs_threshold(yr_threshold = yr_threshold)

yr_threshold<-40
p40<-call_spat_syn_vs_distance_plot_Nyrs_threshold(yr_threshold = yr_threshold)

p40$tab
# A tibble: 10 × 2
# dc                    n
# <fct>             <int>
#   1 (0,25]               36
# 2 (25,75]              53
# 3 (75,150]             59
# 4 (150,250]            59
# 5 (250,400]            68
# 6 (400,700]            74
# 7 (700,1e+03]          71
# 8 (1e+03,1.5e+03]      72
# 9 (1.5e+03,2.5e+03]    60
# 10 (2.5e+03,6e+03]      35

p36$tab
# A tibble: 10 × 2
# dc                    n
# <fct>             <int>
#   1 (0,25]              100
# 2 (25,75]             128
# 3 (75,150]            139
# 4 (150,250]           137
# 5 (250,400]           152
# 6 (400,700]           154
# 7 (700,1e+03]         146
# 8 (1e+03,1.5e+03]     145
# 9 (1.5e+03,2.5e+03]   126
# 10 (2.5e+03,6e+03]      73

p32$tab
# A tibble: 10 × 2
#dc                    n
#<fct>             <int>
#  1 (0,25]              132
#2 (25,75]             191
#3 (75,150]            200
#4 (150,250]           198
#5 (250,400]           215
# 6 (400,700]           208
# 7 (700,1e+03]         205
# 8 (1e+03,1.5e+03]     202
# 9 (1.5e+03,2.5e+03]   172
# 10 (2.5e+03,6e+03]     120
g1<-p40$gp+ylim(0.25,0.8)
g2<-p36$gp+ylim(0.25,0.8)
g3<-p32$gp+ylim(0.25,0.8)

pdf(file = here::here("RESULTS/spat_syn_vs_distance_plot_min40to32yr.pdf"),
    height = 15,
    width = 13)
grid.arrange(g1,g2,g3, nrow=3)
dev.off()










