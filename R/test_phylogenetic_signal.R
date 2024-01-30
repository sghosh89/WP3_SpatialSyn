library(here)
library(ape)
library(phytools)
library(ggplot2)
stree <- read.nexus(here("DATA/BirdTree/tree-pruner-a1bf37e0-0290-4a07-a041-77d2744a7b5c/output.nex"))
# 1000 tree downloaded

#df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_with_optimal_biovar.csv"))
df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_mass.csv"))
dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
df<-left_join(df,dft,by="ScientificName")

df<-df%>%distinct(newBT,.keep_all = T)# just to make sure

# ok, I want to test if any sig phylogenetic signal is there using Pagel's lambda

# initialize
HWI.signal<-data.frame(lambda=NA*numeric(1000),pval=NA*numeric(1000))

for(i in 1:1000){
  tree<-stree[[i]]
  df2<-df %>% arrange(factor(newBT, levels = tree$tip.label))
  nm<-df2$newBT
  
  x2<-df2$HWI
  
  names(x2)<-nm
  
  res2<-phylosig(tree,x=x2,method="lambda",test=TRUE)
  HWI.signal$lambda[i]<-res2$lambda
  HWI.signal$pval[i]<-res2$P
  
  print(i)
}

#res=list(kipps.signal=kipps.signal,
#         HWI.signal=HWI.signal,
#         fLU_ab.signal=fLU_ab.signal)
saveRDS(HWI.signal,here("RESULTS/phylogenetic_signal.RDS"))

# this shows there is phylogenetic signals in traits but not in fLU_ab
hist(HWI.signal$lambda)
hist(HWI.signal$pval)
max(HWI.signal$pval)
sum(HWI.signal$pval<0.05)

phyloHWI<-HWI.signal
sum(phyloHWI$pval<0.05)
g1<-ggplot(phyloHWI,aes(lambda))+geom_histogram(color="black", fill="orange")+
  theme_bw()+xlab(expression("Phylogenetic signal in HWI,"~lambda))

pdf(here("RESULTS/phylogenetic_signal_in_HWI.pdf"), height=3, width=4)
g1
dev.off()



