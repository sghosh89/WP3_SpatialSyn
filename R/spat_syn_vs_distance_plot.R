# plot spatial syn (CorL and CorU) for species against distance category
#rm(list=ls())
library(here)
library(dplyr)
library(tidyverse)
#library(gridExtra)

df<-read.csv(here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_w_morethan2sites.csv"))

dflong<-c()

for(i in 1:nrow(df)){
  
  givenAOU<-df$AOU[i]
  
  distm<-readRDS(here(paste("RESULTS/AOU_",givenAOU,"/distm_sel.RDS",sep="")))
  tempo<-readRDS(here(paste("RESULTS/AOU_",givenAOU,"/abundance_spatsyn/NonParamStat.RDS",sep="")))
  
  indI<-which(tempo$posnI==1,arr.ind = T)
  indN<-which(tempo$posnN==1,arr.ind = T)
  
  cl<-tempo$Corl
  cl[indI]<-NA
  cl[indN]<-NA
  
  cu<-tempo$Coru
  cu[indI]<-NA
  cu[indN]<-NA
  
  #plot(distm,cu)
  
  dtab<-data.frame(d=as.vector(distm),cl=as.vector(cl),cu=as.vector(cu))
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

dflong$dc<-cut(dflong$d, breaks=c(0,250,500,750,1000,1500,3000,4000,8000))

dfs<-dflong%>%group_by(dc)%>%summarise(cl=mean(cl),cu=mean(cu))%>%ungroup()

dg1<-dfs%>%dplyr::select(dc,value=cl)%>%mutate(type="Corl")
dg2<-dfs%>%dplyr::select(dc,value=cu)%>%mutate(type="Coru")
dg<-rbind(dg1,dg2)

g1<-ggplot(dg, aes(x=dc, y=value, col=type)) +
  geom_point()+
  #geom_line()+
  xlab("Pairwise-distance category, Km")+ylab("Spatial synchrony")+
  theme_bw()+theme(legend.position = c(0.7,0.8))

print(g1)

















