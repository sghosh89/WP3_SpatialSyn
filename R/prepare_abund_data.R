# we are preparing data in required format for the copula analysis here 
#===============================================================================
# for a given species ID (i.e. AOU), it will prepare the abundance data in a 
# year by site dataframe format (raw data and detrended data)
# filter applied:
#   - choose only those sites where the species is recorded for atleast 20 out of 23 years
#===============================================================================
#rm(list=ls())
library(here)
library(dplyr)
library(tidyverse)
library(maps)
library(pracma)

prepare_abund_data<-function(givenAOU){
  #chosen_rad<-c(0,400) # radius of interaction, 400 km chosen
  
  resloc<-here(paste("RESULTS/AOU_", givenAOU,sep=""))
  if(!dir.exists(resloc)){
    dir.create(resloc)
  }
  
  xroutes<-read.csv(here("DATA/for_BBS/wrangled_data/uRID_lonlat_stratum.csv"))
  abund_array<-readRDS(here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_array.RDS"))
  #dim(abund_array) # year by site by species
  
  distm<-readRDS(here("DATA/for_BBS/wrangled_data/pairwise_distance_km.RDS"))
  #dim(distm)# site by site
  sitename<-rownames(distm)
  
  spAOU<-unlist(dimnames(abund_array)[3])
  idsp<-which(spAOU%in%givenAOU)
  
  abund_mat<-abund_array[,,idsp]
  
  # now only select sites which are within the interaction radius
  #badid_distm_sel<-which(distm_sel<chosen_rad[1] | distm_sel>chosen_rad[2],arr.ind = T)
  #distm_sel_good<-distm_sel
  #distm_sel_good[badid_distm_sel]<-NA
  #saveRDS(distm_sel_good,paste(resloc,"/distm_sel_good_",chosen_rad[1],"km_",chosen_rad[2],"km.RDS",sep=""))
  
  #----------------------------------------------------------------
  
  # now save abundance matrix in your desired format for the input in tail-analysis
  mat<-(abund_mat) # we will compute tail analysis between sites, i.e., between two columns
  
  # but exclude the sites where the bird species were observed atleast for 20 years
  id<-apply(X=mat, MARGIN=2, FUN = function(x){sum(x==0)})
  badid<-which(id>3)# out of 23 years I only allow 3 absent years
  
  mat_good<-mat[,-badid]
  
  goodsite<-ncol(mat_good)
  
  if(length(badid)>1225){# (=1227-2) means at least 2 sites needed for spatial syn.
    nsites<-NA # no good sites found
  }else{
    #----------- plot a map of birds' present sites (present 20 years at least) --------------
    df_lonlat<-xroutes%>%filter(uRID%in%colnames(mat_good))
    
    wd<-map_data("world")
    wd<-wd%>%filter(region%in%c("USA","Canada"))%>%filter(long<0)
    g1<-ggplot()+coord_fixed()+xlab("")+ylab("")
    g1<-g1+geom_polygon(data=wd, aes(x=long, y=lat, group=group), colour="gray91", fill="gray91")
    g1<-g1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                 panel.background=element_rect(fill="white", colour="white"), axis.line=element_line(colour="white"),
                 legend.position="none",axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
    g1<-g1+geom_point(data=df_lonlat,aes(y=Latitude,x=Longitude),col="black",alpha=0.3,cex=0.3)+
      ggtitle(paste("#sites= ",nrow(df_lonlat),", species' AOU= ",givenAOU,sep=""))+ 
      theme(plot.title = element_text(size = 6, hjust=0.5, vjust=-6))
    g1
    ggsave(paste(resloc,"/routes_on_map_selectedsites.pdf",sep =""),
           width = 10, height = 5, units = "cm")
    #----------------------------------------------------------------
    # save pairwise distance for selected sites
    idsite<-which(sitename%in%colnames(mat_good))
    distm_sel<-distm[idsite,idsite] # distance between present sites
    hist(distm_sel)
    
    saveRDS(distm_sel,paste(resloc,"/distm_sel.RDS",sep=""))
    #-------------------------------------------------------------
    
    mat<-mat_good
    dmat<-pracma::detrend(mat)# detrended matrix: removes linear trend from each column
    
    d_allsite<-vector(mode="list",length=ncol(mat))
    names(d_allsite)<-colnames(mat)
    d_allsite_detrend<-d_allsite
    
    for(i in 1:length(d_allsite)){
      d_allsite[[i]]<-data.frame(Year=rownames(mat),Dat=mat[,i])
      d_allsite_detrend[[i]]<-data.frame(Year=rownames(dmat),Dat=dmat[,i])
    }
    
    saveRDS(mat,paste(resloc,"/year_by_site_abundancedata.RDS",sep=""))
    saveRDS(dmat,paste(resloc,"/year_by_site_abundancedata_detrended.RDS",sep=""))
    
    saveRDS(d_allsite,paste(resloc,"/data_selectedsitelist.RDS",sep=""))
    saveRDS(d_allsite_detrend,paste(resloc,"/detrended_data_selectedsitelist.RDS",sep=""))
    
    #----------------------------------------------------------------------------
    nsites<-length(d_allsite_detrend)
  }
  return(nsites)
}


#=============================================================
# for all the 652 species filtered, get the abundance data
species_absentinfo<-read.csv(here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_absentinfo.csv"))
species_absentinfo$totsites<-1227
#givenAOU<-"5460"
species_absentinfo$ngoodsites<-NA


for(i in 1:nrow(species_absentinfo)){
  givenAOU<-species_absentinfo$AOU[i]
  res<-prepare_abund_data(givenAOU = givenAOU)
  species_absentinfo$ngoodsites[i]<-res
  print(i)
}

write.csv(species_absentinfo,here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_nsites_info.csv"), row.names = F)

# ok, so how many sp are there with good sites?
df<-na.omit(species_absentinfo)
hist(df$ngoodsites,50)

write.csv(df,here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_w_morethan2sites.csv"), row.names = F)
# so 373 species were filtered here based on condition that 
# they were sampled at least at two sites and at min. of 20 years

#df2<-df%>%filter(ngoodsites>250)

