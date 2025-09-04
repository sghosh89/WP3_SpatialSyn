# we are preparing data in required format for the copula analysis here 
#===============================================================================
# for a given species ID (i.e. AOU), it will prepare the abundance data in a 
# year by site dataframe format (raw data and detrended data)
# filter applied:
#   - choose only those sites where the species is recorded for atleast a threshold (vary: 32, 36, 40) out of 41 years
#===============================================================================
#rm(list=ls())
library(here)
library(dplyr)
library(tidyverse)
library(maps)
library(pracma)

prepare_abund_data_Nyrs_threshold<-function(givenAOU, yr_threshold=40){
 
  if(yr_threshold==40){
    resloc<-here(paste("RESULTS/AOU_", givenAOU,sep=""))
    if(!dir.exists(resloc)){
      dir.create(resloc)
    }
    abund_array<-readRDS(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_array.RDS"))
    distm<-readRDS(here("DATA/for_BBS/wrangled_data/pairwise_distance_km_min40yrs.RDS"))
    
  }else{
    resloc2<-here(paste("RESULTS/yr_threshold_", yr_threshold,sep=""))
    if(!dir.exists(resloc2)){
      dir.create(resloc2)
    }
    resloc<-here(paste("RESULTS/yr_threshold_", yr_threshold,"/AOU_", givenAOU,sep=""))
    if(!dir.exists(resloc)){
      dir.create(resloc)
    }
    abund_array<-readRDS(here(paste("DATA/for_BBS/wrangled_data/data1979to2019_abundance_array_min",yr_threshold,"yr.RDS",sep="")))
    distm<-readRDS(here(paste("DATA/for_BBS/wrangled_data/pairwise_distance_km_min",yr_threshold,"yrs.RDS",sep="")))
  }
  
  xroutes<-read.csv(here("DATA/for_BBS/wrangled_data/uRID_lonlat_stratum.csv"))
  
  sitename<-rownames(distm)
  
  spAOU<-unlist(dimnames(abund_array)[3])
  idsp<-which(spAOU%in%givenAOU)
  
  abund_mat<-abund_array[,,idsp]
  
  #----------------------------------------------------------------
  
  # now save abundance matrix in your desired format for the input in tail-analysis
  mat<-(abund_mat) # we will compute tail analysis between sites, i.e., between two columns
  
  # but exclude the sites where the bird species were not observed atleast for 40 years
  id<-apply(X=mat, MARGIN=2, FUN = function(x){sum(x==0)})
  
  # Total number of years in the data (number of rows)
  n_years <- nrow(mat)
  
  # Remove sites where the number of observed years is less than yr_threshold
  badid <- which((n_years - id) < yr_threshold)
  
  mat_good<-mat[,-badid]
  
  goodsite<-ncol(mat_good)
  
  restsites<-dim(abund_array)[2] - 2
  
  if(length(badid)>restsites){# means at least 2 sites needed for spatial syn.
    nsites<-NA # no good sites found
  }else{
    #----------- plot a map of birds' present sites (present 40 years at least) --------------
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
# now call the function:
# yr_threshold<-32; vary it: 32, 36, 40

call_prepare_abund_data_Nyrs_threshold<-function(yr_threshold=40){
  if(yr_threshold==40){
    species_absentinfo<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_absentinfo.csv"))
    abund_array<-readRDS(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_array.RDS"))
    species_absentinfo$totsites<-dim(abund_array)[2]
  }else{
    species_absentinfo<-read.csv(here(paste("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_absentinfo_min",yr_threshold,"yr.csv",sep="")))
    abund_array<-readRDS(here(paste("DATA/for_BBS/wrangled_data/data1979to2019_abundance_array_min",yr_threshold,"yr.RDS",sep="")))
    species_absentinfo$totsites<-dim(abund_array)[2]
  }
  species_absentinfo$ngoodsites<-NA
  
  for(i in 1:nrow(species_absentinfo)){
    givenAOU<-species_absentinfo$AOU[i]
    res<-prepare_abund_data_Nyrs_threshold(givenAOU = givenAOU, yr_threshold=yr_threshold)
    species_absentinfo$ngoodsites[i]<-res
    print(i)
  }
  
  if(yr_threshold==40){
    write.csv(species_absentinfo,here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_nsites_info.csv"), row.names = F)
    
    # ok, so how many sp are there with good sites?
    df<-na.omit(species_absentinfo)
    hist(df$ngoodsites,50)
    as.data.frame(table(df$ngoodsites))# total 173 sp with good sites, out of them 35sp have only 2 good sites
    
    write.csv(df,here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"), row.names = F)
    # so 173 species were filtered here based on condition that 
    # they were sampled at least at two sites and at min. of 40 years
  }else{
    write.csv(species_absentinfo,here(paste("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_nsites_info_min",yr_threshold,"yr.csv", sep="")), row.names = F)
    df<-na.omit(species_absentinfo)
    hist(df$ngoodsites,50)
    as.data.frame(table(df$ngoodsites))# total 173 sp with good sites, out of them 35sp have only 2 good sites
    
    write.csv(df,here(paste("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites_min",yr_threshold,"yr.csv", sep="")), row.names = F)
  }
}

call_prepare_abund_data_Nyrs_threshold(yr_threshold=32)
call_prepare_abund_data_Nyrs_threshold(yr_threshold=36)
call_prepare_abund_data_Nyrs_threshold(yr_threshold = 40)
#df2<-df%>%filter(ngoodsites>250)

