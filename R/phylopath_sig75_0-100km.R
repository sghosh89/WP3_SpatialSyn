rm(list=ls())
library(phylopath);library(RColorBrewer);library(ggtree);library(tidytree)
library(here); library(ape); library(castor); library(phytools)
library(tidyverse); library(gridExtra)
set.seed(seed=123)

if(!dir.exists(here("RESULTS/model_phylopath_sig75_upto100Km"))){
  dir.create(here("RESULTS/model_phylopath_sig75_upto100Km"))
}

if(!dir.exists(here("RESULTS/model_phylopath_sig75_upto100Km/phylopath_LT"))){
  dir.create(here("RESULTS/model_phylopath_sig75_upto100Km/phylopath_LT"))
}

if(!dir.exists(here("RESULTS/model_phylopath_sig75_upto100Km/phylopath_UT"))){
  dir.create(here("RESULTS/model_phylopath_sig75_upto100Km/phylopath_UT"))
}
# ============= read data ==============

df1<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin.csv"))
df1$newBT<-gsub(" ", "_", df1$BirdTreeName)
df1<-df1%>%dplyr::select(AOU,newBT,ScientificName,BirdTreeName)

#-----------
dat<-read.csv(here("RESULTS/distance_sigtaildep_abund_tasmax.csv"))

# split into two category: <=100 km and 100-250 km
df100<-dat%>%filter(dist.KM<=100)
df<-data.frame(AOU=unique(df100$AOU))

# values: total sig within 0-100 Km distance
df$fab.sig<-NA
df$ftasmax.sig<-NA

for(i in 1:nrow(df)){
  
  givenAOU<-df$AOU[i]
  tempo<-df100%>%filter(AOU%in%givenAOU)
  df$fab.sig[i]<-sum(tempo$corlmcoru_ab)/sum(abs(tempo$corlmcoru_ab))
  
  tempoclim<-tempo%>%filter(sig75tasmax==1)
  df$ftasmax.sig[i]<-sum(tempoclim$corlmcoru_ab)/sum(abs(tempoclim$corlmcoru_ab))
}

df$tail<-ifelse(df$fab.sig>0,"LT","UT")
df$tail<-as.factor(df$tail)
df<-left_join(df,df1,by="AOU")
#---------

dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_mass.csv"))
dft<-dft%>%dplyr::select(ScientificName,HWI=meanHWI)
df<-left_join(df,dft,by="ScientificName")

nm<-df
nm<-nm%>%distinct(BirdTreeName)
write.table(nm,here("DATA/BirdTree/unique_speciesnameBirdTree_0_100km_nbin4_tailsig75.txt"),quote=F,col.names =F,row.names=F)

#------------------------------
# remove the duplicated entries from df$new_BT column

df<-df%>%distinct(newBT,.keep_all = T)# just to make sure
df<-df%>%dplyr::select(AOU,fab.sig,ftasmax.sig,HWI,
                       tail,newBT)
df$Species<-df$newBT

sigT<-read.nexus(here("DATA/BirdTree/sig75_0_100km_tree-pruner-9740a40e-4ec7-4417-9940-2461be7e4f2f/output.nex"))
ct<-consensus(sigT, p = 0.5, check.labels = TRUE, rooted = TRUE)# no edge length
ct3<-consensus.edges(sigT,method="least.squares")# with edge length
saveRDS(ct3,here("RESULTS/model_phylopath_sig75_upto100Km/consensus_tree_with_edgelength.RDS"))

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

pdf(here("RESULTS/model_phylopath_sig75_upto100Km/species_phylogeny_0_100km_fLU_ab.pdf"), width = 9, height = 13) # Open a new pdf file
g2 # Write the grid.arrange in the file
dev.off()

pdf(here("RESULTS/model_phylopath_sig75_upto100Km/species_phylogeny_0_100km_HWI.pdf"), width = 9, height = 13) # Open a new pdf file
g3 # Write the grid.arrange in the file
dev.off()

#============== model ===========
modelsHWI_Tonly<-define_model_set(
  model = c(fab.sig~ ftasmax.sig+HWI)
)
gmodels<-plot_model_set(modelsHWI_Tonly, edge_width = 0.5)
ggsave(here("RESULTS/model_phylopath_sig75_upto100Km/modelsHWI_Tonly.pdf"), width=6,height=3)

rownames(df)<-df$Species

#df$ftasmax.sig[which(is.na(df$ftasmax.sig))]<-0

# considering no group
modres_HWI_T_only<- phylo_path(modelsHWI_Tonly, data = df, 
                               tree = ct3, 
                               model = 'lambda')
# 20 sp dropped, 24 sp left

(modsum<-summary(modres_HWI_T_only)) # marginal p-value
gp3<-plot(modsum)+theme_classic()
(best_model_T <- best(modres_HWI_T_only, boot=1000))
gp1<-plot(best_model_T, curvature=0.1, edge_width = 3)
gp2<-coef_plot(best_model_T)+ggplot2::theme_bw()

pdf(here("RESULTS/model_phylopath_sig75_upto100Km/nogroup_model_est.pdf"), height=4, width=10)
grid.arrange(gp1, gp2, ncol=2)
dev.off()

#===================== with groups ===============
dfUT<-df%>%filter(tail=="UT")
dfLT<-df%>%filter(tail=="LT")


# for UT: 11 sp. are dropped from the analysis as they had no sig tail dep synchrony in climate
# 12 data points
sink(here("RESULTS/model_phylopath_sig75_upto100Km/phylopath_UT/modres_HWI_T_only_summary.txt"),
     append=TRUE, split=TRUE)
modres_HWI_T_only<- phylo_path(modelsHWI_Tonly, data = dfUT, 
                               tree = ct3, 
                               model = 'lambda')

(modsum<-summary(modres_HWI_T_only))
gp3<-plot(modsum)+theme_classic()
(best_model_T_UT <- best(modres_HWI_T_only, boot=1000))
gp1<-plot(best_model_T_UT, curvature=0.1, edge_width = 3)
gp2<-coef_plot(best_model_T_UT)+ggplot2::theme_bw()

pdf(here("RESULTS/model_phylopath_sig75_upto100Km/phylopath_UT/model_est.pdf"), height=4, width=10)
grid.arrange(gp1, gp2, ncol=2)
dev.off()

sink()


# for LT: 9 sp. are dropped from the analysis as they had no sig tail dep synchrony in climate
#12 data points
sink(here("RESULTS/model_phylopath_sig75_upto100Km/phylopath_LT/modres_HWI_T_only_summary.txt"),
     append=TRUE, split=TRUE)
modres_HWI_T_only<- phylo_path(modelsHWI_Tonly, data = dfLT, 
                               tree = ct3, 
                               model = 'lambda')

(modsum<-summary(modres_HWI_T_only)) # bad p-values; model fit not well
gp3<-plot(modsum)+theme_classic()
(best_model_T_UT <- best(modres_HWI_T_only, boot=1000))
gp1<-plot(best_model_T_UT, curvature=0.1, edge_width = 3)
gp2<-coef_plot(best_model_T_UT)+ggplot2::theme_bw()

pdf(here("RESULTS/model_phylopath_sig75_upto100Km/phylopath_LT/model_est.pdf"), height=4, width=10)
grid.arrange(gp1, gp2, ncol=2)
dev.off()
sink()





