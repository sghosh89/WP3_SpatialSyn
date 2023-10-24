rm(list=ls())
library(phylopath);library(RColorBrewer);library(ggtree);library(tidytree)
library(here); library(ape); library(castor); library(phytools)
library(tidyverse); library(gridExtra)
set.seed(seed=123)

if(!dir.exists(here("RESULTS/model_phylopath_simple"))){
  dir.create(here("RESULTS/model_phylopath_simple"))
}

if(!dir.exists(here("RESULTS/model_phylopath_simple/phylopath_LT"))){
  dir.create(here("RESULTS/model_phylopath_simple/phylopath_LT"))
}

if(!dir.exists(here("RESULTS/model_phylopath_simple/phylopath_UT"))){
  dir.create(here("RESULTS/model_phylopath_simple/phylopath_UT"))
}
# help: http://blog.phytools.org/2016/03/method-to-compute-consensus-edge.html
wholeT<-read.nexus(here("DATA/BirdTree/whole_tree-pruner-bfb47e7d-3253-4f9e-a5a7-ec93ff54c372/output.nex"))

#=================== first check about your tree description ==========
tree_property<-data.frame(i=1:1000,
                          rooted=NA*numeric(1000),
                          binary=NA*numeric(1000),
                          ultramet=NA*numeric(1000),
                          strict_bifur=NA*numeric(1000))
for (i in c(1:1000)){
  tree<-wholeT[[i]]
  tree_property$rooted[i]<-is.rooted(tree)
  tree_property$binary[i]<-is.binary(tree)
  tree_property$ultramet[i]<-is.ultrametric(tree)
  tree_property$strict_bifur[i]<-castor::is_bifurcating(tree)
}

sum(tree_property$rooted)# should be equal to 1000
sum(tree_property$binary)
sum(tree_property$ultramet)
sum(tree_property$strict_bifur)
# ok, so all fine

# Now get the consensus tree
ct<-consensus(wholeT, p = 0.5, check.labels = TRUE, rooted = TRUE)# no edge length
ct3<-consensus.edges(wholeT,method="least.squares")# with edge length
saveRDS(ct3,here("RESULTS/model_phylopath_simple/consensus_tree_with_edgelength.RDS"))

# ============= read data ==============
df<-read.csv(here("DATA/BirdTree/species_0_250km_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)

dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_with_speciestraits_mass.csv"))
dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
df<-left_join(df,dft,by="ScientificName")

#df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_with_optimal_biovar.csv"))
# it already had the trait kipps and HWI

# remove the duplicated entries from df$new_BT column
df<-df%>%distinct(newBT,.keep_all = T)# just to make sure
df<-df%>%dplyr::select(AOU,fLU_ab,fLU_pr,fLU_tas,kipps, HWI,
                       tail,newBT)
df$Species<-df$newBT
df$tail<-as.factor(df$tail)


#==================== genarate some plots: phylogeny tree =========================
dd<-as_tibble(ct3)
dd<-left_join(dd,df,by=c("label"="newBT"))

g2<-ggtree(ct3,layout="circular") %<+% dd +
  geom_tippoint(pch=19, cex=3,aes(col=fLU_ab))+
  scale_color_gradientn(colours=rev(brewer.pal(n=7,"RdBu")))+
  #scale_color_gradient2(low = "blue",
  #                    midpoint = 0,
  #                    mid = "white",
  #                    high = "red",
  #                    space="Lab", name="Spatial \nsynchrony in \nabundance")+
  #geom_text(aes(label=AOU), hjust=1, vjust=-0.4, size=3)+ 
  theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="gray",
                                               hjust=-0.2)
#g2

g3<-ggtree(ct3,layout="circular") %<+% dd +
  geom_tippoint(pch=19, cex=3,aes(col=HWI))+
  scale_color_gradientn(colours=brewer.pal(n=7,"GnBu"))+
  theme(legend.position="bottom")+ geom_tiplab(aes(label=AOU),color="gray",
                                               hjust=-0.2)
#g3

pdf(here("RESULTS/model_phylopath_simple/species_phylogeny_0_250km_fLU_ab.pdf"), width = 9, height = 13) # Open a new pdf file
g2 # Write the grid.arrange in the file
dev.off()

pdf(here("RESULTS/model_phylopath_simple/species_phylogeny_0_250km_HWI.pdf"), width = 9, height = 13) # Open a new pdf file
g3 # Write the grid.arrange in the file
dev.off()
#=============
dfUT<-df%>%filter(tail=="UT")
dfLT<-df%>%filter(tail=="LT")

rownames(dfUT)<-dfUT$Species
rownames(dfLT)<-dfLT$Species

#=================== Now specify model hypo ===============================
# We are choosing migration ability trait: HWI
modelsHWI_Tonly<-define_model_set(
  model = c(fLU_ab~ fLU_tas+HWI)
)
gmodels<-plot_model_set(modelsHWI_Tonly, edge_width = 0.5)
ggsave(here("RESULTS/model_phylopath_simple/modelsHWI_Tonly.pdf"), width=6,height=3)

#============= test with UT group ================
sink(here("RESULTS/model_phylopath_simple/phylopath_UT/modres_HWI_T_only_summary.txt"),
     append=TRUE, split=TRUE)
modres_HWI_T_only<- phylo_path(modelsHWI_Tonly, data = dfUT, 
                               tree = ct3, 
                               model = 'lambda')

(modsum<-summary(modres_HWI_T_only))
gp3<-plot(modsum)+theme_classic()
(best_model_T_UT <- best(modres_HWI_T_only, boot=1000))
gp1<-plot(best_model_T_UT, curvature=0.1, edge_width = 3)
gp2<-coef_plot(best_model_T_UT)+ggplot2::theme_bw()

pdf(here("RESULTS/model_phylopath_simple/phylopath_UT/model_est.pdf"), height=4, width=10)
grid.arrange(gp1, gp2, ncol=2)
dev.off()

sink()

#============= test with UT group ================
sink(here("RESULTS/model_phylopath_simple/phylopath_LT/modres_HWI_T_only_summary.txt"),
     append=TRUE, split=TRUE)
modres_HWI_T_only<- phylo_path(modelsHWI_Tonly, data = dfLT, 
                               tree = ct3, 
                               model = 'lambda')

(modsum<-summary(modres_HWI_T_only))
gp3<-plot(modsum)+theme_classic()
(best_model_T_UT <- best(modres_HWI_T_only, boot=1000))
gp1<-plot(best_model_T_UT, curvature=0.1, edge_width = 3)
gp2<-coef_plot(best_model_T_UT)+ggplot2::theme_bw()

pdf(here("RESULTS/model_phylopath_simple/phylopath_LT/model_est.pdf"), height=4, width=10)
grid.arrange(gp1, gp2, ncol=2)
dev.off()

sink()



