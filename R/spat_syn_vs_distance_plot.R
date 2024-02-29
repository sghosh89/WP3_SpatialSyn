# plot spatial syn (CorL and CorU) for species against distance category
#rm(list=ls())
library(here)
library(dplyr)
library(tidyverse)
#library(gridExtra)

#df<-read.csv(here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_w_morethan2sites.csv"))
df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))

dflong<-c()

for(i in 1:nrow(df)){
  
  givenAOU<-df$AOU[i]
  
  distm<-readRDS(here(paste("RESULTS/AOU_",givenAOU,"/distm_sel.RDS",sep="")))
  tempo<-readRDS(here(paste("RESULTS/AOU_",givenAOU,"/abundance_spatsyn_nbin_4/NonParamStat.RDS",sep="")))
  
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
  #print(i)
}

# now plot
range(dflong$d)

dflong$dc<-cut(dflong$d, breaks=c(0,25,75,150,250,400,700,1000,1500,2500,5000))

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


#df2s<-dff%>%group_by(dc)%>%summarise(mean.spear=mean(mn.spear),
#                                     min.spear=min(mn.spear),
#                                     max.spear=max(mn.spear))%>%ungroup()


#df2s<-df2%>%group_by(dc)%>%summarise(mn.spear=mean(spear),
#                                     sd.spear=sd(spear))%>%ungroup()
library(PupillometryR)
dff$AOU<-as.factor(dff$AOU)

g2<-ggplot(dff, aes(x=dc, y=mn.spear)) + 
  #geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 1) + 
  geom_point(aes(y = mn.spear, color = AOU), 
             position = position_jitter(width = .15), size = 1, alpha = 0.6) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.2)+ 
  theme_bw()+
  theme(
    #panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    panel.background=element_rect(fill="white", colour="white"), 
    legend.position="none",text=element_text(size=20)
    )+
  ylab("Spatial synchrony in abundance\n averaged across site-pairs")+
  xlab("Between-sites pairwise distance category, Km")#+
  #geom_errorbar(aes(ymin=min.spear, 
   #                 ymax=max.spear), width=.1)+theme_bw() 
  
g2

pdf(here("RESULTS/spat_syn_vs_distance_plot.pdf"),height=5,width = 9)
g2
dev.off()



xx<-dff%>%group_by(dc)%>%summarise(n=n_distinct(AOU))
xx












