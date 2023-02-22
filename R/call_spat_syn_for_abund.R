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
  
  compute_spat_syn(givenAOU = givenAOU, nbin=2, 
                   inputresloc=inputresloc, 
                   outputresloc=outputresloc, target="abundance")
  print(i)
}




