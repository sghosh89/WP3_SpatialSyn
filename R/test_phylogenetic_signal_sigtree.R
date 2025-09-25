library(here)
library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
set.seed(seed=123)
stree40<-read.nexus(here("DATA/BirdTree/sig95_0_250km_tree-pruner-6c06110b-3266-48fc-b9bd-13786bc19ec8/output.nex"))
stree32<-read.nexus(here("DATA/BirdTree/sig95_0_250km_min32yr_tree-pruner-8d5f46d0-b4a6-4c08-a1a1-2aa6c9cd15b0/output.nex"))
stree36<-read.nexus(here("DATA/BirdTree/sig95_0_250km_min36yr_tree-pruner-20e11ea7-5cfc-4cbe-b7ee-055e0b9bdff6/output.nex"))

nbin<-4
yr_threshold<-32# vary this as 32, 36, 40

if(yr_threshold==40){
  stree<-stree40
  dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km.csv",sep="")))
  dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits.csv"))
}else if(yr_threshold==36){
  stree<-stree36
  dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km_min36yr.csv",sep="")))
  dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_min36yr.csv"))
}else{
  stree<-stree32
  dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km_min32yr.csv",sep="")))
  dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_min32yr.csv"))
}


df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin_min32yr.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
df<-df%>%dplyr::select(AOU,newBT,ScientificName,BirdTreeName)

dfsig<-left_join(dfsig,df,by="AOU")

dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
dfsig<-left_join(dfsig,dft,by="ScientificName")

# initialize
df<-dfsig
HWI.signal<-data.frame(lambda=NA*numeric(1000),pval=NA*numeric(1000))
fab.signal<-abs.tot.td.ab.signal<-HWI.signal

for(i in 1:1000){
  tree<-stree[[i]]
  df2<-df %>% arrange(factor(newBT, levels = tree$tip.label))
  nm<-df2$newBT
  
  #----------- for HWI.signal --------------
  x2<-df2$HWI
  names(x2)<-nm
  
  res2<-phylosig(tree,x=x2,method="lambda",test=TRUE)
  HWI.signal$lambda[i]<-res2$lambda
  HWI.signal$pval[i]<-res2$P
  
  #----------- for fab5.signal --------------
  x2<-df2$fab.sig
  names(x2)<-nm
  
  res2<-phylosig(tree,x=x2,method="lambda",test=TRUE)
  fab.signal$lambda[i]<-res2$lambda
  fab.signal$pval[i]<-res2$P
  
  #----------- for abs.tot.td.ab.signal --------------
  x2<-df2$abs.tot.td.ab.sig
  names(x2)<-nm
  
  res2<-phylosig(tree,x=x2,method="lambda",test=TRUE)
  abs.tot.td.ab.signal$lambda[i]<-res2$lambda
  abs.tot.td.ab.signal$pval[i]<-res2$P
  
  print(i)
}

res=list(HWI.signal=HWI.signal,
         fab.signal=fab.signal,
         abs.tot.td.ab.signal=abs.tot.td.ab.signal)

if(yr_threshold==40){
  saveRDS(res,here("RESULTS/phylogenetic_signal_sig95.RDS"))
}else if(yr_threshold==36){
  saveRDS(res,here("RESULTS/phylogenetic_signal_sig95_min36yr.RDS"))
}else{
  saveRDS(res,here("RESULTS/phylogenetic_signal_sig95_min32yr.RDS"))
}

#===================== for 32 years =========================================

res<-readRDS(here("RESULTS/phylogenetic_signal_sig95_min32yr.RDS"))
phyloHWI<-res$HWI.signal#HWI.signal
sum(phyloHWI$pval<0.05)
g1<-ggplot(phyloHWI,aes(lambda))+geom_histogram(color="black", fill="orange")+
  theme_bw()+xlab(expression("Phylogenetic signal in HWI,"~lambda))
g1

mean(phyloHWI$lambda)
mean(res$HWI.signal$lambda)# ~0.99
mean(res$fab.signal$lambda) #~0.0
mean(res$abs.tot.td.ab.signal$lambda) #~0

#===================== for 36 years =========================================
# this shows there is phylogenetic signals in traits but not in fLU_ab

res<-readRDS(here("RESULTS/phylogenetic_signal_sig95_min36yr.RDS"))
phyloHWI<-res$HWI.signal#HWI.signal
sum(phyloHWI$pval<0.05)
g1<-ggplot(phyloHWI,aes(lambda))+geom_histogram(color="black", fill="orange")+
  theme_bw()+xlab(expression("Phylogenetic signal in HWI,"~lambda))
g1

mean(phyloHWI$lambda)
mean(res$HWI.signal$lambda)# ~0.99
mean(res$fab.signal$lambda) #~0.1
mean(res$abs.tot.td.ab.signal$lambda) #~0

#===================== for 40 years =========================================
# this shows there is phylogenetic signals in traits but not in fLU_ab

res<-readRDS(here("RESULTS/phylogenetic_signal_sig95.RDS"))
phyloHWI<-res$HWI.signal#HWI.signal
sum(phyloHWI$pval<0.05)
g1<-ggplot(phyloHWI,aes(lambda))+geom_histogram(color="black", fill="orange")+
  theme_bw()+xlab(expression("Phylogenetic signal in HWI,"~lambda))
g1

mean(phyloHWI$lambda)#0.96
mean(res$HWI.signal$lambda)# ~0.96
mean(res$fab.signal$lambda) #~0.0
mean(res$abs.tot.td.ab.signal$lambda) #~0


