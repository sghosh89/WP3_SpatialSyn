library(here)
library(ape)
stree <- read.nexus(here("DATA/BirdTree/tree-pruner-f1d9b817-3739-4e7d-bbaf-1227c85c4a2c/output.nex"))
# 1000 tree downloaded

# my data table with traits
df<-read.csv(here("DATA/BirdTree/species_0_250km_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_with_speciestraits_mass.csv"))
dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
df<-left_join(df,dft,by="ScientificName")
# remove the duplicated entries from df$new_BT column
df<-df%>%distinct(newBT,.keep_all = T)
df<-df%>%dplyr::select(AOU,fLU_ab,kipps, HWI,tail,newBT)

# ok, I want to test if any sig phylogenetic signal is there using Pagel's lambda
library(phytools)

# initialize
kipps.signal<-data.frame(lambda=NA*numeric(1000),pval=NA*numeric(1000))
HWI.signal<-kipps.signal
fLU_ab.signal<-kipps.signal
  
  
for(i in 1:1000){
  tree<-stree[[i]]
  df2<-df %>% arrange(factor(newBT, levels = tree$tip.label))
  nm<-df2$newBT
  x1<-df2$kipps
  x2<-df2$HWI
  x3<-df2$fLU_ab
  names(x1)<-names(x2)<-names(x3)<-nm
  
  res1<-phylosig(tree,x=x1,method="lambda",test=TRUE)
  kipps.signal$lambda[i]<-res1$lambda
  kipps.signal$pval[i]<-res1$P
  
  res2<-phylosig(tree,x=x2,method="lambda",test=TRUE)
  HWI.signal$lambda[i]<-res2$lambda
  HWI.signal$pval[i]<-res2$P
  
  res3<-phylosig(tree,x=x3,method="lambda",test=TRUE)
  fLU_ab.signal$lambda[i]<-res3$lambda
  fLU_ab.signal$pval[i]<-res3$P
  
  print(i)
}

res=list(kipps.signal=kipps.signal,
         HWI.signal=HWI.signal,
         fLU_ab.signal=fLU_ab.signal)

# this shows there is phylogenetic signals in traits but not in fLU_ab
hist(res$HWI.signal$lambda)
hist(res$HWI.signal$pval)

hist(res$fLU_ab.signal$lambda)
hist(res$fLU_ab.signal$pval)

saveRDS(res,here("RESULTS/phylogenetic_signal.RDS"))






