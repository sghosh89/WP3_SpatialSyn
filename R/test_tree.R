#===================================================================
# We will test some hypotheses related to trait-difference 
# among two groups showing LT and UT spatial synchrony
#===================================================================
#rm(list=ls())
library(readxl)
library(tidyverse)
library(dplyr)
library(here)
library(ggpubr)
library(gridExtra)

df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km.csv"))
df_spmeta<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
df_spmeta<-df_spmeta%>%dplyr::select(AOU, ORDER, Family, Genus, Species, English_Common_Name, ScientificName)

df<-left_join(df,df_spmeta, by="AOU") # don't rename this dataframe
df$BirdTreeName<-NA

#write.csv(df,here("DATA/BirdTree/species_0_250km.csv"),row.names = F)

BT<-read.csv(here("DATA/BirdTree/BLIOCPhyloMasterTax (1).csv"))
idmatch<-which(df$ScientificName%in%BT$Scientific)
df$BirdTreeName[idmatch]<-df$ScientificName[idmatch]

# we will now fill the non-matched name from AVONET Suppmat file/ searching synonyms

id<-which(is.na(df$BirdTreeName)) #70 species
dfnonmatched<-df[id,]
dfnonmatched<-dfnonmatched%>%dplyr::select(ScientificName, BirdTreeName)

Avotalk<-read_excel(here("DATA/AVONET/AVONET Supplementary dataset 1.xlsx"),sheet=11)

dfnonmatched<-left_join(dfnonmatched,Avotalk,by=c("ScientificName"="Species1"))
dfnonmatched$BirdTreeName<-coalesce(dfnonmatched$Species3,dfnonmatched$BirdTreeName)
dfnonmatched<-dfnonmatched%>%dplyr::select(ScientificName, BirdTreeName)

df<-left_join(df,dfnonmatched,by="ScientificName")
df$BirdTreeName<-coalesce(df$BirdTreeName.x,df$BirdTreeName.y)
df<-df%>%dplyr::select(-BirdTreeName.x, -BirdTreeName.y)

write.csv(df,here("DATA/BirdTree/species_0_250km_tobefilled.csv"),row.names = F)

# Now we fill manually the above file column BirdTreeName and saved as
# "DATA/BirdTree/species_0_250km_filledin.csv"


# play with sample data

library(ape)
library(ggtree)
library(tidytree)

stree <- read.nexus(here("DATA/BirdTree/tree-pruner-f1d9b817-3739-4e7d-bbaf-1227c85c4a2c/output.nex"))

ctree<-consensus(stree, p = 1, check.labels = TRUE, rooted = TRUE)

df<-read.csv(here("DATA/BirdTree/species_0_250km_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
df<-df%>%dplyr::select(AOU,fLU_ab,tail,newBT)

dd<-as_tibble(ctree)
dd<-left_join(dd,df,by=c("label"="newBT"))

#tree2<-as.treedata(dd)
#ggtree(ctree, layout="circular")+geom_tiplab()
g1<-ggtree(ctree,layout="circular") %<+% dd +
  geom_tippoint(pch=19, aes(col=tail))+theme(legend.position="bottom")

g2<-ggtree(ctree,layout="circular") %<+% dd +
  geom_tippoint(pch=19, aes(col=fLU_ab))+
  scale_color_gradient2(low = "blue",
                      midpoint = 0,
                      mid = "white",
                      high = "red",
                      space="Lab", name="Spatial \nsynchrony in \nabundance")+
  #geom_text(aes(label=AOU), hjust=1, vjust=-0.4, size=3)+ 
  theme(legend.position="bottom")
#g2

pdf(here("RESULTS/species_phylogeny_0_250km.pdf"), width = 10, height = 6) # Open a new pdf file
grid.arrange(g1, g2, nrow=1) # Write the grid.arrange in the file
dev.off()


