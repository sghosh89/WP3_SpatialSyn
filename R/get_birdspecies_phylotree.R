#rm(list=ls())
library(readxl)
library(tidyverse)
library(dplyr)
library(here)
library(ggpubr)
library(gridExtra)

#library(ape)
#library(ggtree)

df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km.csv"))
df_spmeta<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
df_spmeta<-df_spmeta%>%dplyr::select(AOU, ORDER, Family, Genus, Species, English_Common_Name, ScientificName)

df<-left_join(df,df_spmeta, by="AOU") # don't rename this dataframe
df$BirdTreeName<-NA

#write.csv(df,here("DATA/BirdTree/species_0_250km.csv"),row.names = F)

BT<-read.csv(here("DATA/BirdTree/BLIOCPhyloMasterTax.csv"))
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

df<-read.csv(here("DATA/BirdTree/species_0_250km_filledin.csv"))
nm<-df
nm<-nm%>%distinct(BirdTreeName)
write.table(nm,here("DATA/BirdTree/unique_speciesnameBirdTree_0_250km.txt"),quote=F,col.names =F,row.names=F)
# 254 unique sp name in BirdTree

################ get species phylogeny for LT and UT separate group ###########
df<-read.csv(here("DATA/BirdTree/species_0_250km_filledin.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)

dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_with_speciestraits_mass.csv"))
dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
df<-left_join(df,dft,by="ScientificName")

# remove the duplicated entries from df$new_BT column
df<-df%>%distinct(newBT,.keep_all = T)

dfl<-df%>%filter(tail=="LT")# 124 unique sp.
btl<-dfl%>%dplyr::select(BirdTreeName)
write.table(btl,here("DATA/BirdTree/unique_speciesnameBirdTree_0_250km_LTabund.txt"),quote=F,col.names =F,row.names=F)

dfu<-df%>%filter(tail=="UT") #130 unique species
btu<-dfu%>%dplyr::select(BirdTreeName)
write.table(btu,here("DATA/BirdTree/unique_speciesnameBirdTree_0_250km_UTabund.txt"),quote=F,col.names =F,row.names=F)

