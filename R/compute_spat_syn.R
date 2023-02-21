# This script will do the copula analysis, will generate necessary plots
# for a given metapopulation or env. variables across the sites

library(here)
library(dplyr)
library(tidyverse)

source(here("R/NonParamStat.R"))
source(here("R/NonParamStat_matrixplot.R"))

# Input:
# mat: a matrix or dataframe where each target species time series along each column, 
#  rows have name for sampling years
#resloc: path to save results
#nbin: 2 (default) to measure tail-dep.

# Output:
# a list and several plots to be saved in resloc path

compute_spat_syn<-function(givenAOU, nbin=2){
  
  inputresloc<-here(paste("RESULTS/AOU_", givenAOU,sep=""))
  outputresloc<-here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn",sep=""))
  
  if(!dir.exists(outputresloc)){
    dir.create(outputresloc)
  }
  
  distm<-readRDS(here(paste(inputresloc,"/distm_sel.RDS",sep="")))
  
  d_allsite_detrend<-readRDS(here(paste(inputresloc,"/detrended_data_selectedsitelist.RDS",sep="")))
  
  # just to check
  #distm<-distm[1:50,1:50]
  #d_allsite_detrend<-d_allsite_detrend[1:50]
  #nbin<-2
  
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

#givenAOU<-"5460"

df<-read.csv(here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_w_morethan2sites.csv"))
#for(i in 1:3){
for(i in 1:nrow(df)){
  givenAOU<-df$AOU[i]
  compute_spat_syn(givenAOU = givenAOU, nbin=2)
  print(i)
}




