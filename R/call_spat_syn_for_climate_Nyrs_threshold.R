# This script will compute spat syn among climate timeseries from different sites
library(here)
source(here("R/compute_spat_syn.R"))
#============================================
# Now call for detrended climate data 

call_spat_syn_for_climate_Nyrs_threshold<-function(climvar, nbin=4, yr_threshold=40){
  
  if(yr_threshold==40){
    df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))
  }else{
    df<-read.csv(here(paste("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites_min",yr_threshold,"yr.csv",sep="")))
  }
  
  for(i in 1:nrow(df)){
    
    givenAOU<-df$AOU[i]
    
    if(yr_threshold==40){
      inputresloc<-here(paste("RESULTS/AOU_",givenAOU,sep=""))
      
      if(climvar=="tas_avgAprtoAug"){
        d_allsite_detrend<-readRDS(here(paste(inputresloc,"/tas_detrended_data_avgAprtoAug_selectedsitelist.RDS",sep="")))
        #}else if(climvar=="tas_avgMaytoJuly"){
        #  d_allsite_detrend<-readRDS(here(paste(inputresloc,"/tas_detrended_data_avgMaytoJuly_selectedsitelist.RDS",sep="")))
        #}else if(climvar=="tas_avgAprtoJuly"){
        #  d_allsite_detrend<-readRDS(here(paste(inputresloc,"/tas_detrended_data_avgAprtoJuly_selectedsitelist.RDS",sep="")))
      }else if(climvar=="pr_avgAprtoAug"){
        d_allsite_detrend<-readRDS(here(paste(inputresloc,"/pr_detrended_data_avgAprtoAug_selectedsitelist.RDS",sep="")))
      }else{
        d_allsite_detrend<-readRDS(here(paste(inputresloc,"/",climvar,"_detrended_data_selectedsitelist.RDS",sep="")))
      }
      
    }else{
      inputresloc<-here(paste("RESULTS/yr_threshold_",yr_threshold,"/AOU_",givenAOU,sep=""))
      
      if(climvar=="tas_avgAprtoAug"){
        d_allsite_detrend<-readRDS(here(paste(inputresloc,"/tas_detrended_data_avgAprtoAug_selectedsitelist.RDS",sep="")))
        #}else if(climvar=="tas_avgMaytoJuly"){
        #  d_allsite_detrend<-readRDS(here(paste(inputresloc,"/tas_detrended_data_avgMaytoJuly_selectedsitelist.RDS",sep="")))
        #}else if(climvar=="tas_avgAprtoJuly"){
        #  d_allsite_detrend<-readRDS(here(paste(inputresloc,"/tas_detrended_data_avgAprtoJuly_selectedsitelist.RDS",sep="")))
      }else if(climvar=="pr_avgAprtoAug"){
        d_allsite_detrend<-readRDS(here(paste(inputresloc,"/pr_detrended_data_avgAprtoAug_selectedsitelist.RDS",sep="")))
      }else{
        d_allsite_detrend<-readRDS(here(paste(inputresloc,"/",climvar,"_detrended_data_selectedsitelist.RDS",sep="")))
      }
      
    }
    
    distm<-readRDS(here(paste(inputresloc,"/distm_sel.RDS",sep="")))
    
    #===================================
    #nbin<-4
    
    if(yr_threshold==40){
      outputresloc<-here(paste("RESULTS/AOU_", givenAOU,"/",climvar,"_spatsyn_nbin_",nbin,sep=""))
    }else{
      outputresloc<-here(paste("RESULTS/yr_threshold_",yr_threshold,"/AOU_", givenAOU,"/",climvar,"_spatsyn_nbin_",nbin,sep=""))
    }
    
    
    if(!dir.exists(outputresloc)){
      dir.create(outputresloc)
    }
    
    compute_spat_syn(nbin=nbin, 
                     inputresloc=inputresloc, 
                     outputresloc=outputresloc, distm=distm, d_allsite=d_allsite_detrend)
    
    print(i)
    
  }
  cat("-------done---------\n")
}

# call the function
climvar<-"pr"
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 32)
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 36)
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 40)

climvar<-"pr_avgAprtoAug"
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 32)
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 36)
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 40)

climvar<-"tas"
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 32)
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 36)
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 40)

climvar<-"tas_avgAprtoAug"
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 32)
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 36)
call_spat_syn_for_climate_Nyrs_threshold(climvar=climvar,nbin=4, yr_threshold = 40)

