rm(list=ls())
library(phylopath);library(RColorBrewer);library(ggtree);library(tidytree)
library(here); library(ape); library(castor); library(phytools)
library(tidyverse); library(gridExtra)
set.seed(seed=123)

if(!dir.exists(here("RESULTS/model_phylopath_sig75_100-250km"))){
  dir.create(here("RESULTS/model_phylopath_sig75_100-250km"))
}
# ============= read data ==============

df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
df<-df%>%dplyr::select(AOU,newBT,ScientificName,BirdTreeName)

#-----------
nbin<-4
dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail75sig_summary_100-250Km.csv",sep="")))
dfsig<-left_join(dfsig,df,by="AOU")
#---------

dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_mass.csv"))
dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
dfsig<-left_join(dfsig,dft,by="ScientificName")

nm<-dfsig
nm<-nm%>%distinct(BirdTreeName)
write.table(nm,here("DATA/BirdTree/unique_speciesnameBirdTree_100_250km_nbin4_tailsig75.txt"),quote=F,col.names =F,row.names=F)
# 50 sp.
#------------------------------
# remove the duplicated entries from df$new_BT column
df<-dfsig
df<-df%>%distinct(newBT,.keep_all = T)# just to make sure
df<-df%>%dplyr::select(AOU,
                       abs.tot.td.ab.sig,
                       tot.td.ab.sig,
                       fab.sig,
                       ftas.sig,ftas5.sig,
                       abs.tot.td.tas5.sig,
                       tot.td.tas5.sig,
                       ftas4.sig,ftas3.sig,
                       fpr.sig,
                       fpr5.sig,
                       abs.tot.td.pr5.sig,
                       tot.td.pr5.sig,
                       HWI,
                       tail75,newBT)
df$Species<-df$newBT
df$tail75<-as.factor(df$tail75)

sigT<-read.nexus(here("DATA/BirdTree/sig75_100_250km_tree-pruner-b22209ef-1dc3-4b39-a54a-2ec23a92f78d/output.nex"))
ct<-consensus(sigT, p = 0.5, check.labels = TRUE, rooted = TRUE)# no edge length
ct3<-consensus.edges(sigT,method="least.squares")# with edge length
saveRDS(ct3,here("RESULTS/model_phylopath_sig75_100-250Km/consensus_tree_with_edgelength.RDS"))

dd<-as_tibble(ct3)
dd<-left_join(dd,df,by=c("label"="newBT"))
g2<-ggtree(ct3,layout="circular") %<+% dd +
  geom_tippoint(pch=19, cex=3,aes(col=fab.sig))+
  scale_color_gradientn(colours=rev(brewer.pal(n=5,"RdBu")))+
  #scale_color_gradient2(low = "blue",
  #                    midpoint = 0,
  #                    mid = "white",
  #                    high = "red",
  #                    space="Lab", name="Spatial \nsynchrony in \nabundance")+
  #geom_text(aes(label=AOU), hjust=1, vjust=-0.4, size=3)+ 
  theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="gray",
                                               hjust=-0.2)
g2

g3<-ggtree(ct3,layout="circular") %<+% dd +
  geom_tippoint(pch=19, cex=3,aes(col=HWI))+
  scale_color_gradientn(colours=brewer.pal(n=5,"GnBu"))+
  theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="gray",
                                               hjust=-0.2)
g3

pdf(here("RESULTS/model_phylopath_sig75_100-250Km/species_phylogeny_0_100km_fLU_ab.pdf"), width = 9, height = 13) # Open a new pdf file
g2 # Write the grid.arrange in the file
dev.off()

pdf(here("RESULTS/model_phylopath_sig75_100-250Km/species_phylogeny_0_100km_HWI.pdf"), width = 9, height = 13) # Open a new pdf file
g3 # Write the grid.arrange in the file
dev.off()

#============== model ===========
call_phylopath_sig75_100_250km<-function(model, df){
  
  #-------------------------------
  model_tas5_abs_td<-define_model_set(
    model = c(abs.tot.td.ab.sig~ abs.tot.td.tas5.sig+HWI)
  )
  
  model_pr5_abs_td<-define_model_set(
    model = c(abs.tot.td.ab.sig~ abs.tot.td.pr5.sig+HWI)
  )
  
  model_tas5_tot_td<-define_model_set(
    model = c(tot.td.ab.sig~ tot.td.tas5.sig+HWI)
  )
  
  model_pr5_tot_td<-define_model_set(
    model = c(tot.td.ab.sig~ tot.td.pr5.sig+HWI)
  )
  
  model_tas<-define_model_set(
    model = c(fab.sig~ ftas.sig+HWI)
  )
  
  model_tas5<-define_model_set(
    model = c(fab.sig~ ftas5.sig+HWI)
  )
  
  model_tas4<-define_model_set(
    model = c(fab.sig~ ftas4.sig+HWI)
  )
  
  model_tas3<-define_model_set(
    model = c(fab.sig~ ftas3.sig+HWI)
  )
  
  model_pr<-define_model_set(
    model = c(fab.sig~ fpr.sig+HWI)
  )
  
  model_pr5<-define_model_set(
    model = c(fab.sig~ fpr5.sig+HWI)
  )
  
  if(model=="tas"){
    modelsHWI_Tonly<-model_tas
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_tas")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else if(model=="tas5"){
    modelsHWI_Tonly<-model_tas5
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_tas5")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else if(model=="tas4"){
    modelsHWI_Tonly<-model_tas4
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_tas4")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else if(model=="tas3"){
    modelsHWI_Tonly<-model_tas3
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_tas3")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else if(model=="pr"){
    modelsHWI_Tonly<-model_pr
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_pr")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else if(model=="pr5"){
    modelsHWI_Tonly<-model_pr5
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_pr5")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else if(model=="tas5_abs_td"){
    modelsHWI_Tonly<-model_tas5_abs_td
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_tas5_abs_td")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else if(model=="pr5_abs_td"){
    modelsHWI_Tonly<-model_pr5_abs_td
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_pr5_abs_td")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else if(model=="tas5_tot_td"){
    modelsHWI_Tonly<-model_tas5_tot_td
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_tas5_tot_td")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else if(model=="pr5_tot_td"){
    modelsHWI_Tonly<-model_pr5_tot_td
    myresloc<-here("RESULTS/model_phylopath_sig75_100-250km/model_pr5_tot_td")
    if(!dir.exists(myresloc)){dir.create(myresloc)}
  }else{
    print("define model")
  }
  #-------------------------
  gmodels<-plot_model_set(modelsHWI_Tonly, edge_width = 0.5)
  ggsave(here(paste(myresloc,"modelsHWI_Tonly.pdf",sep="/")), width=6,height=3)
  
  rownames(df)<-df$Species
  
  #------------- no group --------------
  sink(here(paste(myresloc,"modres_HWI_T_only_summary.txt",sep="/")),
       append=TRUE, split=TRUE)
  
  modres_HWI_T_only<- phylo_path(modelsHWI_Tonly, data = df, 
                                 tree = ct3, 
                                 model = 'lambda')
  
  modsum<-summary(modres_HWI_T_only)
  print(modsum)
  gp3<-plot(modsum)+theme_classic()
  best_model_T <- best(modres_HWI_T_only, boot=1000)
  print(best_model_T)
  gp1<-plot(best_model_T, curvature=0.1, edge_width = 3)
  gp2<-coef_plot(best_model_T)+ggplot2::theme_bw()
  
  sink()
  
  pdf(here(paste(myresloc,"nogroup_model_est.pdf",sep="/")), height=4, width=10)
  grid.arrange(gp1, gp2, ncol=2)
  dev.off()
  
  #------------- with groups --------------
  
  dfUT<-df%>%filter(tail75=="UT")# 
  dfLT<-df%>%filter(tail75=="LT")# 
  
  myreslocUT<-here(paste(myresloc,"UT",sep="/"))
  myreslocLT<-here(paste(myresloc,"LT",sep="/"))
  
  if(!dir.exists(myreslocUT)){dir.create(myreslocUT)}
  if(!dir.exists(myreslocLT)){dir.create(myreslocLT)}
  
  #------------- with UT --------------
  sink(here(paste(myreslocUT,"modres_HWI_T_only_summary.txt",sep="/")),
       append=TRUE, split=TRUE)
  modres_HWI_UT_only<- phylo_path(modelsHWI_Tonly, data = dfUT, 
                                  tree = ct3, 
                                  model = 'lambda')
  
  modsumUT<-summary(modres_HWI_UT_only)
  print(modsumUT)
  gp3<-plot(modsumUT)+theme_classic()
  best_model_T_UT <- best(modres_HWI_UT_only, boot=1000)
  print(best_model_T_UT)
  gp1<-plot(best_model_T_UT, curvature=0.1, edge_width = 3)
  gp2<-coef_plot(best_model_T_UT)+ggplot2::theme_bw()
  
  pdf(here(paste(myreslocUT,"model_est.pdf",sep="/")), height=4, width=10)
  grid.arrange(gp1, gp2, ncol=2)
  dev.off()
  
  sink()
  #------------- with LT --------------
  sink(here(paste(myreslocLT,"modres_HWI_T_only_summary.txt",sep="/")),
       append=TRUE, split=TRUE)
  modres_HWI_LT_only<- phylo_path(modelsHWI_Tonly, data = dfLT, 
                                  tree = ct3, 
                                  model = 'lambda')
  
  modsumLT<-summary(modres_HWI_LT_only)
  print(modsumLT)
  gp3<-plot(modsumLT)+theme_classic()
  best_model_T_LT <- best(modres_HWI_LT_only, boot=1000)
  print(best_model_T_LT)
  gp1<-plot(best_model_T_LT, curvature=0.1, edge_width = 3)
  gp2<-coef_plot(best_model_T_LT)+ggplot2::theme_bw()
  
  pdf(here(paste(myreslocLT,"model_est.pdf",sep="/")), height=4, width=10)
  grid.arrange(gp1, gp2, ncol=2)
  dev.off()
  sink()
}
#=========================

call_phylopath_sig75_100_250km(model="tas", df=df)
call_phylopath_sig75_100_250km(model="tas5", df=df)
call_phylopath_sig75_100_250km(model="tas3", df=df)
call_phylopath_sig75_100_250km(model="tas4", df=df)

call_phylopath_sig75_100_250km(model="pr", df=df)
call_phylopath_sig75_100_250km(model="pr5", df=df)

call_phylopath_sig75_100_250km(model="tas5_abs_td", df=df)
call_phylopath_sig75_100_250km(model="tas5_tot_td", df=df)

call_phylopath_sig75_100_250km(model="pr5_abs_td", df=df)
call_phylopath_sig75_100_250km(model="pr5_tot_td", df=df)
