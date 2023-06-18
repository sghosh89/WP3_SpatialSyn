# We want to see which area on the map shows more LT across species and which area shows more UT
rm(list=ls())
library(here)
library(tidyverse)
library(maps)
library(gridExtra)
df_spmeta<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
df_spmeta<-df_spmeta%>%dplyr::select(AOU,ngoodsites,"ScientificName")

df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km.csv"))
df<-left_join(df,df_spmeta,by="AOU")
# select LT in abundance, and atleast >80% for species
dfLT<-df%>%filter(tail=="LT")%>% filter(abs(fLU_ab)>0.9)

gp_list_LT <- vector(mode='list', length=nrow(dfLT))
# now plot on map for those 27 species the sites they are sampled

xroutes<-read.csv(here("DATA/for_BBS/wrangled_data/uRID_lonlat_stratum.csv"))

for(i in 1:nrow(dfLT)){
  aou<-dfLT$AOU[i]
  resloc<-here(paste("RESULTS/AOU_",aou,sep=""))
  readfile<-readRDS(paste(resloc,"/year_by_site_abundancedata_detrended.RDS",sep=""))
  sites<-colnames(readfile)
  mylonlat<-xroutes%>%filter(uRID %in% sites)
  
  wd<-map_data("world")
  wd<-wd%>%filter(region%in%c("USA","Canada"))%>%filter(long<0)
  g1<-ggplot()+coord_fixed()+xlab("")+ylab("")
  g1<-g1+geom_polygon(data=wd, aes(x=long, y=lat, group=group), colour="white", fill="white")
  g1<-g1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
               panel.background=element_rect(fill="gray", colour="white"), 
               axis.line=element_line(colour="white"),
               legend.position="none",
               axis.ticks=element_blank(), 
               axis.text.x=element_blank(), 
               axis.text.y=element_blank())
  g1<-g1+geom_point(data=mylonlat,aes(y=Latitude,x=Longitude),col="red",alpha=0.3,cex=1)+
    #ggtitle(paste("#sites= ",nrow(mylonlat),", species' AOU= ",aou,sep=""))+ 
    ggtitle(paste("AOU= ",aou,sep=""))+
    theme(plot.title = element_text(size = 12, hjust=0.5, vjust=0))
  gp_list_LT[[i]]<-g1
  
}


pdf(here("RESULTS/allmapsLT.pdf"),width=10,height=7)
grid.arrange(gp_list_LT[[1]], gp_list_LT[[2]], gp_list_LT[[3]],
             gp_list_LT[[4]], gp_list_LT[[5]], gp_list_LT[[6]],
             gp_list_LT[[7]], gp_list_LT[[8]], gp_list_LT[[9]],
             gp_list_LT[[10]], gp_list_LT[[11]], gp_list_LT[[12]],
             gp_list_LT[[13]], gp_list_LT[[14]], gp_list_LT[[15]],
             gp_list_LT[[16]], gp_list_LT[[17]], gp_list_LT[[18]],
             gp_list_LT[[19]], gp_list_LT[[20]], gp_list_LT[[21]],
             gp_list_LT[[22]], gp_list_LT[[23]], gp_list_LT[[24]],
             gp_list_LT[[25]], gp_list_LT[[26]], gp_list_LT[[27]],
             ncol=6)
dev.off()

# AOU = 4881, 5275, 7315 was excluded as not identified unambiguously

# did not plot for UT, only plotted for LT

