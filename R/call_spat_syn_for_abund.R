# This script will compute spat syn among metapopulation's abundance
library(here)
source(here("R/compute_spat_syn.R"))
#============================================
# Now call for species' abundance data 

df<-read.csv(here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_w_morethan2sites.csv"))

for(i in 1:nrow(df)){
  
  givenAOU<-df$AOU[i]
  
  inputresloc<-here(paste("RESULTS/AOU_", givenAOU,sep=""))
  outputresloc<-here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn",sep=""))
  
  if(!dir.exists(outputresloc)){
    dir.create(outputresloc)
  }
  
  distm<-readRDS(here(paste(inputresloc,"/distm_sel.RDS",sep="")))
  d_allsite_detrend<-readRDS(here(paste(inputresloc,"/detrended_data_selectedsitelist.RDS",sep="")))
  
  compute_spat_syn(nbin=2, 
                   inputresloc=inputresloc, 
                   outputresloc=outputresloc, distm=distm, d_allsite=d_allsite_detrend)
  print(i)
}




