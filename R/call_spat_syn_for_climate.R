# This script will compute spat syn among climate timeseries from different sites
library(here)
source(here("R/compute_spat_syn.R"))
#============================================
# Now call for species' abundance data 


call_spat_syn_for_climate<-function(climvar){
  
  df<-read.csv(here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_w_morethan2sites.csv"))
  for(i in 1:nrow(df)){
    
    givenAOU<-df$AOU[i]
    
    inputresloc<-here(paste("RESULTS/AOU_", givenAOU,sep=""))
    outputresloc<-here(paste("RESULTS/AOU_", givenAOU,"/",climvar,"_spatsyn",sep=""))
    
    if(!dir.exists(outputresloc)){
      dir.create(outputresloc)
    }
    
    distm<-readRDS(here(paste(inputresloc,"/distm_sel.RDS",sep="")))
    d_allsite_detrend<-readRDS(here(paste(inputresloc,"/",climvar,"_detrended_data_selectedsitelist.RDS",sep="")))
    
    compute_spat_syn(nbin=2, 
                     inputresloc=inputresloc, 
                     outputresloc=outputresloc, distm=distm, d_allsite=d_allsite_detrend)
    print(i)
  }
  cat("-------done---------\n")
}

# call the function
climvar<-"pr"
call_spat_syn_for_climate(climvar=climvar)

climvar<-"tas"
call_spat_syn_for_climate(climvar=climvar)

climvar<-"tasmax"
call_spat_syn_for_climate(climvar=climvar)

climvar<-"tasmin"
call_spat_syn_for_climate(climvar=climvar)

