rm(list=ls())
library(phylolm);library(RColorBrewer);library(ggtree);library(tidytree)
library(here); library(ape); library(castor); library(phytools); library(caper)
library(tidyverse); library(gridExtra)
set.seed(seed=123)

# ============= read data ==============

df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
df<-df%>%dplyr::select(AOU,newBT,ScientificName,BirdTreeName)

nbin<-4
dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km.csv",sep="")))
dfsig<-left_join(dfsig,df,by="AOU")


dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_mass.csv"))
dft<-dft%>%dplyr::select(ScientificName,HWI=meanHWI, Diet.5Cat)
dfsig<-left_join(dfsig,dft,by="ScientificName")
dfsig$Diet.5Cat<-as.factor(dfsig$Diet.5Cat)

# remove the duplicated entries from df$new_BT column
df<-dfsig
df<-df%>%distinct(newBT,.keep_all = T)# just to make sure
df<-df%>%dplyr::select(AOU,
                       abs.tot.td.ab.sig,
                       fab.sig,
                       Diet.5Cat,
                       tail95,newBT)
df$Species<-df$newBT
df$tail95<-as.factor(df$tail95)

diet.gr<-df$Diet.5Cat
names(diet.gr)<-df$Species
sigT<-read.nexus(here("DATA/BirdTree/sig95_0_250km_tree-pruner-6c06110b-3266-48fc-b9bd-13786bc19ec8/output.nex"))
ct<-consensus(sigT, p = 0.5, check.labels = TRUE, rooted = TRUE)# no edge length
ct3<-consensus.edges(sigT,method="least.squares")# with edge length
fab.sig<-df$fab.sig
names(fab.sig)<-df$Species
abs.tot.td.ab.sig<-df$abs.tot.td.ab.sig
names(abs.tot.td.ab.sig)<-df$Species
phyl_phytools<-phylANOVA(ct3,diet.gr,fab.sig, nsim=1000,posthoc=TRUE, p.adj=c("holm"))
print(phyl_phytools)

#=====================================
# read and print coefficients
#xx<-readRDS(here("RESULTS/model_phylolm_sig75_0-250km/model_tas5/model_est_phylolm_boot.RDS"))
#y<-summary(xx)
#round(y$coefficients,4)
#round(y$bootmean,4)
#round(y$bootconfint95,4)
#round(y$r.squared, 4); round(y$adj.r.squared, 4)
#y






