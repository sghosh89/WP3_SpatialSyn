# This script will do the copula analysis, will generate necessary plots
# for a given metapopulation or env. variables across the sites

library(here)
library(dplyr)
library(tidyverse)

source(here("R/NonParamStat.R"))
source(here("R/NonParamStat_matrixplot.R"))

# Input:
#givenAOU = unique species ID for BBS data
#nbin: 2 (default) to measure tail-dep.
# inputresloc: path for the location from where 
#               it should read the input list (e.g.,species abundance list of each site 
#                 when target=="abundance")
# outputresloc: path for the location where it should save the result
# target: either "abundance", or "temperature" 

# Output:
# save results for spatial syn at either tail

compute_spat_syn<-function(givenAOU, nbin=2, inputresloc, outputresloc, target="abundance"){
  
  if(target=="abundance"){
    
    distm<-readRDS(here(paste(inputresloc,"/distm_sel.RDS",sep="")))
    d_allsite_detrend<-readRDS(here(paste(inputresloc,"/detrended_data_selectedsitelist.RDS",sep="")))
    
    resloc<-paste(outputresloc,"/",sep="")
    
    z<-multcall(d_allsite = d_allsite_detrend, resloc = resloc, nbin=nbin)
    
    saveRDS(z,paste(resloc,"NonParamStat.RDS",sep=""))
    NonParamStat_matrixplot(data=z,
                            resloc=resloc,
                            tl.cex=1.2,cl.cex=2,line=1)
    
    spear<-z$spear # this is all positive correlation (sometimes -ve values but they were indep), as invert the -ve corr copula
    posnN<-z$posnN
    posnI<-z$posnI
    
    # plot spat syn (spearman correlation) vs distm 
    corval<-z$corval
    corval[posnI]<-NA
    #distm[posnI]<-NA # you don't need it for the below plot
    
    png(paste(resloc,file="Corval_vs_distance.png",sep=''), width = 2000,height = 2000,  res = 300)
    plot(as.vector(distm),as.vector(corval),type="p",xlab="pairwise distance, Km",
         ylab="pairwise Spearman correlation", pch=16, col=rgb(0,0,0,0.2)) # this would be a decreasing relationship: spat syn should decrease as distance increase
    dev.off()
    
  }

  #
  #if(target=="temperature"){
    # fill this in to compute spat syn among temperature timeseries
  #}
  
}


