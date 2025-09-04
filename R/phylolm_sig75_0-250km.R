rm(list=ls())
library(phylolm);library(RColorBrewer);library(ggtree);library(tidytree)
library(here); library(ape); library(castor); library(phytools); library(caper)
library(tidyverse); library(gridExtra)
set.seed(seed=123)

if(!dir.exists(here("RESULTS/model_phylolm_sig75_0-250km"))){
  dir.create(here("RESULTS/model_phylolm_sig75_0-250km"))
}

# ============= read data ==============

df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
df<-df%>%dplyr::select(AOU,newBT,ScientificName,BirdTreeName)

nbin<-4
dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail75sig_summary_0-250Km.csv",sep="")))
dfsig<-left_join(dfsig,df,by="AOU")


dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits.csv"))
dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
dfsig<-left_join(dfsig,dft,by="ScientificName")

nm<-dfsig
nm<-nm%>%distinct(BirdTreeName)
write.table(nm,here("DATA/BirdTree/unique_speciesnameBirdTree_0_250km_nbin4_tailsig75.txt"),quote=F,col.names =F,row.names=F)
# Now we will download 1000 trees with the above 59 species names


# remove the duplicated entries from df$new_BT column
df<-dfsig
df<-df%>%distinct(newBT,.keep_all = T)# just to make sure
df<-df%>%dplyr::select(AOU,
                       abs.tot.td.ab.sig,
                       fab.sig,
                       ftas.sig,
                       ftas5.sig,
                       abs.tot.td.tas5.sig,
                       fpr.sig,
                       fpr5.sig,
                       abs.tot.td.pr5.sig,
                       HWI,
                       tail75,newBT)
df$Species<-df$newBT
df$tail75<-as.factor(df$tail75)


sigT<-read.nexus(here("DATA/BirdTree/sig75_0_250km_tree-pruner-dd4d5870-1c32-4bfb-9253-d3af12edd234/output.nex"))
ct<-consensus(sigT, p = 0.5, check.labels = TRUE, rooted = TRUE)# no edge length
ct3<-consensus.edges(sigT,method="least.squares")# with edge length
saveRDS(ct3,here("RESULTS/model_phylolm_sig75_0-250km/consensus_tree_with_edgelength.RDS"))


dd<-as_tibble(ct3)
dd<-left_join(dd,df,by=c("label"="newBT"))

g1<-ggtree(ct3,layout="circular") %<+% dd +
  geom_tippoint(pch=19, cex=6,aes(col=abs.tot.td.ab.sig))+
  scale_color_gradientn(colours=brewer.pal(n=5,"PuRd"))+
  theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="black",
                                               hjust=-0.3)
g1

g2<-ggtree(ct3,layout="circular") %<+% dd +
  geom_tippoint(pch=19, cex=6,aes(col=fab.sig))+
  scale_color_gradientn(colours=rev(brewer.pal(n=5,"RdBu")))+
  theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="black",
                                               hjust=-0.3)
g2

g3<-ggtree(ct3,layout="circular") %<+% dd +
  geom_tippoint(pch=19, cex=6,aes(col=HWI))+
  scale_color_gradientn(colours=brewer.pal(n=5,"GnBu"))+
  theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="black",
                                               hjust=-0.3)
#g3

pdf(here("RESULTS/model_phylolm_sig75_0-250km/species_phylogeny_0_250km_absolute_ab.pdf"), width = 9, height = 9) # Open a new pdf file
g1 # Write the grid.arrange in the file
dev.off()

pdf(here("RESULTS/model_phylolm_sig75_0-250km/species_phylogeny_0_250km_fLU_ab.pdf"), width = 9, height = 9) # Open a new pdf file
g2 # Write the grid.arrange in the file
dev.off()

pdf(here("RESULTS/model_phylolm_sig75_0-250km/species_phylogeny_0_250km_HWI.pdf"), width = 9, height = 9) # Open a new pdf file
g3 # Write the grid.arrange in the file
dev.off()

#============== model ===========
call_phylolm_sig75_0_250km<-function(model, df, ct3){
  
  set.seed(seed=123)
  
  rownames(df)<-df$Species
  ct3null<-ct3
  ct3null$node.label<-NULL
  
  if(model=="tas"){
    fit <- phylolm(fab.sig~ ftas.sig+HWI,data=df,
                   phy=ct3,model="lambda")
    fitboot <- phylolm(fab.sig~ ftas.sig+HWI,data=df,
                       phy=ct3,model="lambda")
    myresloc<-here("RESULTS/model_phylolm_sig75_0-250km/model_tas")
    
    fit_pgls<-pgls(fab.sig~ ftas.sig+HWI,
                   comparative.data(ct3null,df,"Species"), 
                   lambda="ML", bounds = list(lambda=c(1e-06,1),
                                              kappa = c(1e-06,3), 
                                              delta = c(1e-06,3)))
    
    if(!dir.exists(myresloc)){dir.create(myresloc)}
    saveRDS(fit,here(paste(myresloc,"model_est_phylolm.RDS",sep="/")))
    saveRDS(fitboot,here(paste(myresloc,"model_est_phylolm_boot.RDS",sep="/")))
    saveRDS(fit_pgls,here(paste(myresloc,"model_est_pgls.RDS",sep="/")))
    
  }else if(model=="tas5"){
    
    fit <- phylolm(fab.sig~ ftas5.sig+HWI,data=df,
                   phy=ct3,model="lambda")
    fitboot <- phylolm(fab.sig~ ftas5.sig+HWI,data=df,
                       phy=ct3,model="lambda", boot=1000)
    
    myresloc<-here("RESULTS/model_phylolm_sig75_0-250km/model_tas5")
    fit_pgls<-pgls(fab.sig~ ftas5.sig+HWI,
                   comparative.data(ct3null,df,"Species"), 
                   lambda="ML", bounds = list(lambda=c(1e-06,1),
                                              kappa = c(1e-06,3), 
                                              delta = c(1e-06,3)))
    
    if(!dir.exists(myresloc)){dir.create(myresloc)}
    saveRDS(fit,here(paste(myresloc,"model_est_phylolm.RDS",sep="/")))
    saveRDS(fitboot,here(paste(myresloc,"model_est_phylolm_boot.RDS",sep="/")))
    saveRDS(fit_pgls,here(paste(myresloc,"model_est_pgls.RDS",sep="/")))
    
  }else if(model=="pr"){
    
    fit <- phylolm(fab.sig~ fpr.sig+HWI,data=df,
                   phy=ct3,model="lambda")
    fitboot <- phylolm(fab.sig~ fpr.sig+HWI,data=df,
                       phy=ct3,model="lambda", boot=1000)
    
    myresloc<-here("RESULTS/model_phylolm_sig75_0-250km/model_pr")
    fit_pgls<-pgls(fab.sig~ fpr.sig+HWI,
                   comparative.data(ct3null,df,"Species"), 
                   lambda="ML", bounds = list(lambda=c(1e-06,1),
                                              kappa = c(1e-06,3), 
                                              delta = c(1e-06,3)))
    
    if(!dir.exists(myresloc)){dir.create(myresloc)}
    saveRDS(fit,here(paste(myresloc,"model_est_phylolm.RDS",sep="/")))
    saveRDS(fitboot,here(paste(myresloc,"model_est_phylolm_boot.RDS",sep="/")))
    saveRDS(fit_pgls,here(paste(myresloc,"model_est_pgls.RDS",sep="/")))
    
  }else if(model=="pr5"){
    
    fit <- phylolm(fab.sig~ fpr5.sig+HWI,data=df,
                   phy=ct3,model="lambda")
    fitboot <- phylolm(fab.sig~ fpr5.sig+HWI,data=df,
                       phy=ct3,model="lambda", boot=1000)
    myresloc<-here("RESULTS/model_phylolm_sig75_0-250km/model_pr5")
    fit_pgls<-pgls(fab.sig~ fpr5.sig+HWI,
                   comparative.data(ct3null,df,"Species"), 
                   lambda="ML", bounds = list(lambda=c(1e-06,1),
                                              kappa = c(1e-06,3), 
                                              delta = c(1e-06,3)))
    
    if(!dir.exists(myresloc)){dir.create(myresloc)}
    saveRDS(fit,here(paste(myresloc,"model_est_phylolm.RDS",sep="/")))
    saveRDS(fitboot,here(paste(myresloc,"model_est_phylolm_boot.RDS",sep="/")))
    saveRDS(fit_pgls,here(paste(myresloc,"model_est_pgls.RDS",sep="/")))
    
  }else if(model=="tas5_abs_td"){
    
    fit <- phylolm(abs.tot.td.ab.sig~ abs.tot.td.tas5.sig+HWI,data=df,
                   phy=ct3,model="lambda")
    fitboot <- phylolm(abs.tot.td.ab.sig~ abs.tot.td.tas5.sig+HWI,data=df,
                       phy=ct3,model="lambda", boot=1000)
    
    myresloc<-here("RESULTS/model_phylolm_sig75_0-250km/model_tas5_abs_td")
    fit_pgls<-pgls(abs.tot.td.ab.sig~ abs.tot.td.tas5.sig+HWI,
                   comparative.data(ct3null,df,"Species"), 
                   lambda="ML", bounds = list(lambda=c(1e-06,1),
                                              kappa = c(1e-06,3), 
                                              delta = c(1e-06,3)))
    
    if(!dir.exists(myresloc)){dir.create(myresloc)}
    saveRDS(fit,here(paste(myresloc,"model_est_phylolm.RDS",sep="/")))
    saveRDS(fitboot,here(paste(myresloc,"model_est_phylolm_boot.RDS",sep="/")))
    saveRDS(fit_pgls,here(paste(myresloc,"model_est_pgls.RDS",sep="/")))
    
  }else if(model=="pr5_abs_td"){
    
    fit <- phylolm(abs.tot.td.ab.sig~ abs.tot.td.pr5.sig+HWI,data=df,
                   phy=ct3,model="lambda")
    fitboot <- phylolm(abs.tot.td.ab.sig~ abs.tot.td.pr5.sig+HWI,data=df,
                       phy=ct3,model="lambda", boot=1000)
    myresloc<-here("RESULTS/model_phylolm_sig75_0-250km/model_pr5_abs_td")
    fit_pgls<-pgls(abs.tot.td.ab.sig~ abs.tot.td.pr5.sig+HWI,
                   comparative.data(ct3null,df,"Species"), 
                   lambda="ML", bounds = list(lambda=c(1e-06,1),
                                              kappa = c(1e-06,3), 
                                              delta = c(1e-06,3)))
    
    if(!dir.exists(myresloc)){dir.create(myresloc)}
    saveRDS(fit,here(paste(myresloc,"model_est_phylolm.RDS",sep="/")))
    saveRDS(fitboot,here(paste(myresloc,"model_est_phylolm_boot.RDS",sep="/")))
    saveRDS(fit_pgls,here(paste(myresloc,"model_est_pgls.RDS",sep="/")))
    
  }else{
    print("define model")
  }
  #-------------------------
  
  sink(here(paste(myresloc,"modres_HWI_T_only_summary.txt",sep="/")),
       append=TRUE, split=TRUE)
  
  cat("============== results using phylolm function ============\n")
  print(summary(fit))
  
  cat("============== results using phylolm function with bootstrap ============\n")
  print(summary(fitboot))
  
  cat("\n============== results using pgls function ============\n")
  print(summary(fit_pgls))
  sink()
}
#=====================================================
call_phylolm_sig75_0_250km(model="tas5_abs_td", df=df, ct3=ct3)
call_phylolm_sig75_0_250km(model="pr5_abs_td", df=df, ct3=ct3)

call_phylolm_sig75_0_250km(model="tas5", df=df, ct3=ct3)
call_phylolm_sig75_0_250km(model="pr5", df=df, ct3=ct3)

call_phylolm_sig75_0_250km(model="tas", df=df, ct3=ct3)
call_phylolm_sig75_0_250km(model="pr", df=df, ct3=ct3)


