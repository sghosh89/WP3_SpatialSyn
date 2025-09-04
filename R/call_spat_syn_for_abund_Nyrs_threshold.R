# This script will compute spat syn among metapopulation's abundance
library(here)
source(here("R/compute_spat_syn.R"))
#============================================
# Now call for species' abundance data 

call_spat_syn_for_abund_Nyrs_threshold<-function(df,nbin, yr_threshold=40){
  
  for(i in 1:nrow(df)){
    
    givenAOU<-df$AOU[i]
    
    if(yr_threshold==40){
      inputresloc<-here(paste("RESULTS/AOU_", givenAOU,sep=""))
      outputresloc<-here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn_nbin_",nbin,sep=""))
      if(!dir.exists(outputresloc)){
        dir.create(outputresloc)
      }
    }else{
      inputresloc<-here(paste("RESULTS/yr_threshold_",yr_threshold,"/AOU_",givenAOU,sep=""))
      outputresloc<-here(paste("RESULTS/yr_threshold_",yr_threshold,"/AOU_",givenAOU,"/abundance_spatsyn_nbin_",nbin,sep=""))
      if(!dir.exists(outputresloc)){
        dir.create(outputresloc)
      }
    }
    
    distm<-readRDS(here(paste(inputresloc,"/distm_sel.RDS",sep="")))
    d_allsite_detrend<-readRDS(here(paste(inputresloc,"/detrended_data_selectedsitelist.RDS",sep="")))
    
    #===================================
    #nbin<-4
    
    compute_spat_syn(nbin=nbin, # for min 40 data points nbin=4 is okay
                     inputresloc=inputresloc, 
                     outputresloc=outputresloc, distm=distm, d_allsite=d_allsite_detrend)
    
    #=============================
    print(i)
  }
  
}


# ok, now call
df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites_min32yr.csv"))
call_spat_syn_for_abund_Nyrs_threshold(df =df, nbin=2, yr_threshold = 32)
call_spat_syn_for_abund_Nyrs_threshold(df =df, nbin=4, yr_threshold = 32)

df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites_min36yr.csv"))
call_spat_syn_for_abund_Nyrs_threshold(df =df, nbin=2, yr_threshold = 36)
call_spat_syn_for_abund_Nyrs_threshold(df =df, nbin=4, yr_threshold = 36)

df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))
call_spat_syn_for_abund_Nyrs_threshold(df =df, nbin=2, yr_threshold = 40)
call_spat_syn_for_abund_Nyrs_threshold(df =df, nbin=4, yr_threshold = 40)
