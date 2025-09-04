#rm(list=ls())
library(readxl)
library(tidyverse)
library(dplyr)
library(here)
library(ggpubr)
library(gridExtra)


# read whole species
#mytrait<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_min32yr.csv"))

df32<-read.csv(here("RESULTS/abundance_spatsyn_nbin_4_tail75sig_summary_0-250Km_min32yr.csv"))
#df32.1<-read.csv(here("RESULTS/abundance_spatsyn_nbin_4_tail95sig_summary_0-250Km_min32yr.csv"))
#df32.1$AOU%in%df32$AOU
#df<-read.csv(here("RESULTS/abundance_spatsyn_nbin_4_tail75sig_summary_0-250Km.csv"))
#df36$AOU%in%df32$AOU
#df$AOU%in%df32$AOU


df_spmeta<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.csv"))  
df_spmeta<-df_spmeta%>%dplyr::select(AOU, ORDER, Family, Genus, Species, English_Common_Name, ScientificName)

df<-left_join(df32,df_spmeta, by="AOU") # don't rename this dataframe
df$BirdTreeName<-NA

#write.csv(df,here("DATA/BirdTree/species_0_250km.csv"),row.names = F)

BT<-read.csv(here("DATA/BirdTree/BLIOCPhyloMasterTax.csv"))
idmatch<-which(df$ScientificName%in%BT$Scientific)
df$BirdTreeName[idmatch]<-df$ScientificName[idmatch]

# we will now fill the non-matched name from AVONET Suppmat file/ searching synonyms

id<-which(is.na(df$BirdTreeName)) #48 species
dfnonmatched<-df[id,]
dfnonmatched<-dfnonmatched%>%dplyr::select(ScientificName, BirdTreeName)

Avotalk<-read_excel(here("DATA/AVONET/AVONET Supplementary dataset 1.xlsx"),sheet=11)

dfnonmatched<-left_join(dfnonmatched,Avotalk,by=c("ScientificName"="Species1"))
dfnonmatched$BirdTreeName<-coalesce(dfnonmatched$Species3,dfnonmatched$BirdTreeName)
dfnonmatched<-dfnonmatched%>%dplyr::select(ScientificName, BirdTreeName)

df<-left_join(df,dfnonmatched,by="ScientificName")
df$BirdTreeName<-coalesce(df$BirdTreeName.x,df$BirdTreeName.y)
df<-df%>%dplyr::select(-BirdTreeName.x, -BirdTreeName.y)

df$BirdTreeName[which(df$ScientificName=="Colaptes auratus auratus")]<-"Colaptes auratus"

# ok, see here AOU=6520 occurred twice with two different BirdTree names (Dendroica aestiva and Dendroica petechia) for a given 
# Scientific species (Setophaga petechia). I think we should exclude the Dendroica aestiva.

id<-which(df$BirdTreeName=="Dendroica aestiva")
df<-df[-id,]
# now it is 184 species!

df$BirdTreeName[which(df$ScientificName=="Setophaga coronata coronata")]<-"Dendroica coronata" # matched with english common name:Yellow-rumped Warbler
df$BirdTreeName[which(df$ScientificName=="Junco hyemalis hyemalis")]<-"Junco hyemalis" # matched with english common name: dark-eyed junco
df$BirdTreeName[which(df$ScientificName=="Leucophaeus atricilla")]<-"Larus atricilla" # matched with english common name: Laughing Gull

# below subsp. were not found in BirdTREE database, so filling with possible sp.
df$BirdTreeName[which(df$ScientificName=="Junco hyemalis oreganus")]<-"Junco hyemalis"
df$BirdTreeName[which(df$ScientificName=="Colaptes auratus cafer")]<-"Colaptes auratus"
df$BirdTreeName[which(df$ScientificName=="Setophaga coronata audoboni")]<-"Dendroica coronata"

write.csv(df,here("DATA/BirdTree/species_0_250km_nbin_4_filledin_min32yr.csv"),row.names = F)

#df<-read.csv(here("DATA/BirdTree/species_0_250km_filledin.csv"))
nm<-df
nm<-nm%>%distinct(BirdTreeName)
write.table(nm,here("DATA/BirdTree/unique_speciesnameBirdTree_0_250km_nbin_4_min32yr.txt"),quote=F,col.names =F,row.names=F)
# 180 unique sp name in BirdTree, 4 sp. are species = subspecies level

