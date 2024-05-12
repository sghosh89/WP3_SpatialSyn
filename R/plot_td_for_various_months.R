rm(list=ls())
library(tidyverse)
library(here)
library(gridExtra)
#----------------

plot_td_various_months<-function(nbin=4, target_dist_cat, siglevel="75"){
  
  dff<-readRDS(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,
                          "_corlmcoru_sigres_summary_",target_dist_cat[1],"-",
                          target_dist_cat[2],"Km.RDS",sep="")))
  
  if(siglevel=="75"){
    dff<-dff%>%dplyr::select(AOU, Lsig75ab, Usig75ab,
                             # significant total lower and upper tail-dep. 
                             # between distance category for climate
                             # considered exactly thoese site pairs for which
                             # we got sig. tail-dep in abundance
                             # significance level: 75%
                             Lsig75tas, Usig75tas,
                             Lsig75tas5, Usig75tas5,
                             Lsig75pr, Usig75pr,
                             Lsig75pr5, Usig75pr5)
    
    
    dff<-dff%>%filter(Lsig75ab!=0 | Usig75ab!=0) # keep species showing sig. tail-dep spatial synchrony
    
    dff<-dff%>%dplyr::select(-Lsig75ab, -Usig75ab)
    
    dff<-gather(dff, key="td", value="measure", -AOU)
    dd<-dff%>%group_by(AOU, td)%>%summarise(mn.measure=mean(measure))%>%ungroup()
    dd$td<-as.factor(dd$td)
    dd$AOU<-as.character(dd$AOU)
    
    ddpr<-dd%>%filter(td%in%c("Lsig75pr", "Usig75pr",
                              "Lsig75pr5", "Usig75pr5"))
    ddtas<-dd%>%filter(td%in%c("Lsig75tas", "Usig75tas",
                               "Lsig75tas5", "Usig75tas5"))
  }
  
  if(siglevel=="95"){
    dff<-dff%>%dplyr::select(AOU, Lsig95ab, Usig95ab,
                             # significant total lower and upper tail-dep. 
                             # between distance category for climate
                             # considered exactly thoese site pairs for which
                             # we got sig. tail-dep in abundance
                             # significance level: 95%
                             Lsig95tas, Usig95tas,
                             Lsig95tas5, Usig95tas5,
                             Lsig95pr, Usig95pr,
                             Lsig95pr5, Usig95pr5)
    
    
    dff<-dff%>%filter(Lsig95ab!=0 | Usig95ab!=0) # keep species showing sig. tail-dep spatial synchrony
    
    dff<-dff%>%dplyr::select(-Lsig95ab, -Usig95ab)
    
    dff<-gather(dff, key="td", value="measure", -AOU)
    dd<-dff%>%group_by(AOU, td)%>%summarise(mn.measure=mean(measure))%>%ungroup()
    dd$td<-as.factor(dd$td)
    dd$AOU<-as.character(dd$AOU)
    
    ddpr<-dd%>%filter(td%in%c("Lsig95pr", "Usig95pr",
                              "Lsig95pr5", "Usig95pr5"))
    ddtas<-dd%>%filter(td%in%c("Lsig95tas", "Usig95tas",
                               "Lsig95tas5", "Usig95tas5"))
  }
  
  gpr<-ggplot(data=ddpr, aes(x=AOU, y=mn.measure, fill=td), alpha=0.5) +
    geom_bar(position="stack", stat="identity")+
    xlab("Bird species, AOU")+
    ylab(paste("Tail-dependence in precipitation for site-pairs within",
               target_dist_cat[1],"-",target_dist_cat[2],"Km",sep=""))+
    theme_bw()+theme(legend.position = "none")+coord_flip()+
    scale_fill_manual(values = c("#FDDBC7", "#B2182B","#D1E5F0", "#2166AC"))
  
  gtas<-ggplot(data=ddtas, aes(x=AOU, y=mn.measure, fill=td)) +
    geom_bar(position="stack", stat="identity")+
    xlab("")+
    ylab(paste("Tail-dependence in temperature for site-pairs within",
         target_dist_cat[1],"-",target_dist_cat[2],"Km",sep=""))+
    theme_bw()+theme(legend.position = c(0.8, 0.2), 
                     axis.text.y=element_blank(),axis.ticks.y=element_blank())+coord_flip()+
    #scale_fill_brewer(palette="RdBu")
    scale_fill_manual(values = c("#FDDBC7", "#B2182B","#D1E5F0", "#2166AC"))
  
  pdf(here(paste("RESULTS/plot_td_for_various_months_",
                 target_dist_cat[1],"-",target_dist_cat[2],"Km_siglevel_",siglevel,".pdf",sep="")), width=10, height=10)
  grid.arrange(gpr,gtas, ncol=2)
  dev.off()
  
}


target_dist_cat<-c(0,250)
plot_td_various_months(nbin=4, target_dist_cat=target_dist_cat, siglevel = "75")

target_dist_cat<-c(0,250)
plot_td_various_months(nbin=4, target_dist_cat=target_dist_cat, siglevel = "95")

#target_dist_cat<-c(0,100)
#plot_td_various_months(nbin=4, target_dist_cat=target_dist_cat)

#target_dist_cat<-c(100,250)
#plot_td_various_months(nbin=4, target_dist_cat=target_dist_cat)




