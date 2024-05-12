library(here)
library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
set.seed(seed=123)

#stree<-read.nexus(here("DATA/BirdTree/sig75_0_250km_tree-pruner-dd4d5870-1c32-4bfb-9253-d3af12edd234/output.nex"))
stree<-read.nexus(here("DATA/BirdTree/sig95_0_250km_tree-pruner-6c06110b-3266-48fc-b9bd-13786bc19ec8/output.nex"))

df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
df<-df%>%dplyr::select(AOU,newBT,ScientificName,BirdTreeName)

nbin<-4
#dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail75sig_summary_0-250Km.csv",sep="")))
dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km.csv",sep="")))
dfsig<-left_join(dfsig,df,by="AOU")


dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_mass.csv"))
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
#saveRDS(res,here("RESULTS/phylogenetic_signal_sig75.RDS"))
saveRDS(res,here("RESULTS/phylogenetic_signal_sig95.RDS"))

# this shows there is phylogenetic signals in traits but not in fLU_ab

phyloHWI<-res$fab.signal#HWI.signal
sum(phyloHWI$pval<0.05)
g1<-ggplot(phyloHWI,aes(lambda))+geom_histogram(color="black", fill="orange")+
  theme_bw()+xlab(expression("Phylogenetic signal in HWI,"~lambda))
g1

mean(phyloHWI$lambda)

#pdf(here("RESULTS/phylogenetic_signal_sig75_in_HWI.pdf"), height=3, width=4)
#g1
#dev.off()
res<-readRDS(here("RESULTS/phylogenetic_signal_sig95.RDS"))
mean(res$HWI.signal$lambda)# ~0.96
mean(res$fab.signal$lambda) #~0
mean(res$abs.tot.td.ab.signal$lambda) #~0



