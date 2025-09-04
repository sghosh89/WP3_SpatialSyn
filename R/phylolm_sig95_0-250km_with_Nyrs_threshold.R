rm(list=ls())
library(phylolm);library(RColorBrewer);library(ggtree);library(tidytree)
library(here); library(ape); library(castor); library(phytools); library(caper)
library(tidyverse); library(gridExtra)


#yr_threshold<-32

get_df<-function(yr_threshold){
  
  set.seed(seed=123)
  
  if(yr_threshold==40){
    outputresloc<-here("RESULTS/model_phylolm_sig95_0-250km")
    if(!dir.exists(outputresloc)){
      dir.create(outputresloc)
    }
  }else{
    outputresloc<-here(paste0("RESULTS/model_phylolm_sig95_0-250km_min",yr_threshold,"yr"))
    if(!dir.exists(outputresloc)){
      dir.create(outputresloc)
    }
  }
  
  
  # ============= read data ==============
  
  df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin_min32yr.csv"))
  df$newBT<-gsub(" ", "_", df$BirdTreeName)
  df<-df%>%dplyr::select(AOU,newBT,ScientificName,BirdTreeName)
  
  nbin<-4
  
  if(yr_threshold==40){
    dfsig<-read.csv(here(paste0("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km.csv")))
    dft<-read.csv(here(paste0("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_",nbin,"_with_speciestraits.csv")))
  }else{
    dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km_min",yr_threshold,"yr.csv",sep="")))
    dft<-read.csv(here(paste("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_",nbin,"_with_speciestraits_min",yr_threshold,"yr.csv",sep="")))
  }
  
  dfsig<-left_join(dfsig,df,by="AOU")
  
  dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
  dfsig<-left_join(dfsig,dft,by="ScientificName")
  
  nm<-dfsig
  #xx<-dfsig[dfsig$BirdTreeName %in% dfsig$BirdTreeName[duplicated(dfsig$BirdTreeName)], ] # common English name: Pacific Wren appears as Winter Wren in BirdTree
  nm<-nm%>%distinct(BirdTreeName) # 134 unique species in BirdTree
  
  if(yr_threshold==40){
    write.table(nm,here("DATA/BirdTree/unique_speciesnameBirdTree_0_250km_nbin4_tailsig95.txt"),quote=F,col.names =F,row.names=F)
    # Now we will download 1000 trees with the above 35 species names
  }else{
    write.table(nm,here(paste("DATA/BirdTree/unique_speciesnameBirdTree_0_250km_nbin4_tailsig95_min",yr_threshold,"yr.txt",sep="")),quote=F,col.names =F,row.names=F)
    # Now we will download 1000 trees with the above 134 species names (for 32 yr)
  }
  
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
                         tail95,newBT)
  df$Species<-df$newBT
  df$tail95<-as.factor(df$tail95)
  saveRDS(df,here(paste(outputresloc,"/df.RDS",sep="")))
  
  ############################################################
  
  if(yr_threshold==32){
    sigT<-read.nexus(here(paste("DATA/BirdTree/sig95_0_250km_min",yr_threshold,"yr_tree-pruner-8d5f46d0-b4a6-4c08-a1a1-2aa6c9cd15b0/output.nex",sep="")))
  }
  
  if(yr_threshold==36){
    sigT<-read.nexus(here(paste("DATA/BirdTree/sig95_0_250km_min",yr_threshold,"yr_tree-pruner-20e11ea7-5cfc-4cbe-b7ee-055e0b9bdff6/output.nex",sep="")))
  }
  
  if(yr_threshold==40){
    sigT<-read.nexus(here("DATA/BirdTree/sig95_0_250km_tree-pruner-6c06110b-3266-48fc-b9bd-13786bc19ec8/output.nex"))
  }
  
  ct<-consensus(sigT, p = 0.5, check.labels = TRUE, rooted = TRUE)# no edge length
  ct3<-consensus.edges(sigT,method="least.squares")# with edge length
  
  saveRDS(ct3,here(paste(outputresloc,"/consensus_tree_with_edgelength.RDS",sep="")))
  
  dd<-as_tibble(ct3)
  dd<-left_join(dd,df,by=c("label"="newBT"))
  
  g1<-ggtree(ct3,layout="circular") %<+% dd +
    geom_tippoint(pch=19, cex=6,aes(col=abs.tot.td.ab.sig))+
    scale_color_gradientn(colours=brewer.pal(n=5,"PuRd"))+
    theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="black",
                                                 hjust=-0.3)
  g2<-ggtree(ct3,layout="circular") %<+% dd +
    geom_tippoint(pch=19, cex=6,aes(col=fab.sig))+
    scale_color_gradientn(colours=rev(brewer.pal(n=5,"RdBu")))+
    theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="black",
                                                 hjust=-0.3)
  g3<-ggtree(ct3,layout="circular") %<+% dd +
    geom_tippoint(pch=19, cex=6,aes(col=HWI))+
    scale_color_gradientn(colours=brewer.pal(n=5,"GnBu"))+
    theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="black",
                                                 hjust=-0.3)
  
  pdf(here(paste(outputresloc,"/species_phylogeny_0_250km_absolute_ab.pdf",sep="")), width = 9, height = 9) # Open a new pdf file
  print(g1) # Write the grid.arrange in the file
  dev.off()
  
  pdf(here(paste(outputresloc,"/species_phylogeny_0_250km_fLU_ab.pdf",sep="")), width = 9, height = 9) # Open a new pdf file
  print(g2) # Write the grid.arrange in the file
  dev.off()
  
  pdf(here(paste(outputresloc,"/species_phylogeny_0_250km_HWI.pdf",sep="")), width = 9, height = 9) # Open a new pdf file
  print(g3) # Write the grid.arrange in the file
  dev.off()
}

#============== model ===========
call_phylolm_sig95_0_250km_Nyrs_threshold<-function(model, yr_threshold){
  
  set.seed(seed=123)
  
  if(yr_threshold==40){
    outputresloc<-here("RESULTS/model_phylolm_sig95_0-250km")
  }else{
    outputresloc<-here(paste0("RESULTS/model_phylolm_sig95_0-250km_min",yr_threshold,"yr"))
  }
  
  ct3<-readRDS(here(paste(outputresloc,"/consensus_tree_with_edgelength.RDS",sep="")))
  
  ct3null<-ct3
  ct3null$node.label<-NULL
  
  df<-readRDS(here(paste(outputresloc,"/df.RDS",sep="")))
  rownames(df)<-df$Species
  
  if(model=="tas"){
    fit <- phylolm(fab.sig~ ftas.sig+HWI,data=df,
                   phy=ct3,model="lambda")
    fitboot <- phylolm(fab.sig~ ftas.sig+HWI,data=df,
                       phy=ct3,model="lambda")
    myresloc<-here(paste0(outputresloc,"/model_tas"))
    
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
    
    myresloc<-here(paste0(outputresloc,"/model_tas5"))
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
    
    myresloc<-here(paste0(outputresloc,"/model_pr"))
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
    myresloc<-here(paste0(outputresloc,"/model_pr5"))
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
    
    myresloc<-here(paste0(outputresloc,"/model_tas5_abs_td"))
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
    myresloc<-here(paste0(outputresloc,"/model_pr5_abs_td"))
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

get_df(yr_threshold = 32)

call_phylolm_sig95_0_250km_Nyrs_threshold(model="tas5_abs_td", yr_threshold=32)
call_phylolm_sig95_0_250km_Nyrs_threshold(model="pr5_abs_td", yr_threshold=32)

call_phylolm_sig95_0_250km_Nyrs_threshold(model="tas5", yr_threshold=32)
call_phylolm_sig95_0_250km_Nyrs_threshold(model="pr5", yr_threshold=32)

call_phylolm_sig95_0_250km_Nyrs_threshold(model="tas", yr_threshold=32)
call_phylolm_sig95_0_250km_Nyrs_threshold(model="pr", yr_threshold=32)

#=====================================================

get_df(yr_threshold = 36)

call_phylolm_sig95_0_250km_Nyrs_threshold(model="tas5_abs_td", yr_threshold=36)
call_phylolm_sig95_0_250km_Nyrs_threshold(model="pr5_abs_td", yr_threshold=36)

call_phylolm_sig95_0_250km_Nyrs_threshold(model="tas5", yr_threshold=36)
call_phylolm_sig95_0_250km_Nyrs_threshold(model="pr5", yr_threshold=36)

call_phylolm_sig95_0_250km_Nyrs_threshold(model="tas", yr_threshold=36)
call_phylolm_sig95_0_250km_Nyrs_threshold(model="pr", yr_threshold=36)


#=====================================================

get_df(yr_threshold = 40)

call_phylolm_sig95_0_250km_Nyrs_threshold(model="tas5_abs_td", yr_threshold=40)
call_phylolm_sig95_0_250km_Nyrs_threshold(model="pr5_abs_td", yr_threshold=40)

call_phylolm_sig95_0_250km_Nyrs_threshold(model="tas5", yr_threshold=40)
call_phylolm_sig95_0_250km_Nyrs_threshold(model="pr5", yr_threshold=40)

call_phylolm_sig95_0_250km_Nyrs_threshold(model="tas", yr_threshold=40)
call_phylolm_sig95_0_250km_Nyrs_threshold(model="pr", yr_threshold=40)




