#====================================================================
# This file will get the body size related data for the bird species
# considered for spatial synchrony (0-250Km) results, 
#====================================================================

#rm(list=ls())
`%notin%` <- Negate(`%in%`)
library(tidyverse)
library(dplyr)
library(readxl)
library(here)
# first for all bird species

avodat_meta<-read_excel(here("DATA/AVONET/AVONET Supplementary dataset 1.xlsx"),sheet=1)
avodat_raw<-read_excel(here("DATA/AVONET/AVONET Supplementary dataset 1.xlsx"),sheet=6)
#avodat_mass =? separate data sheets - needed to be combined

#=========================================================================================================
df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_min32yr.csv")) # 32 yrs. = yr_threshold
#=========================================================================================================

df_spmeta<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.csv"))
df_spmeta<-df_spmeta%>%dplyr::select(AOU, English_Common_Name, ScientificName)
df<-left_join(df,df_spmeta, by="AOU")

# ok, see the problem
# species names are not same throughout the three sources of avodat_raw
# Species1_BirdLife, Species2_eBird, Species3_BirdTree
#
mydat<-data.frame(spname=df$ScientificName)
# first match with source2
get_dat_ebird<-left_join(mydat,avodat_raw,by=c("spname"="Species2_eBird"))
length(unique(get_dat_ebird$Avibase.ID)) #250

id<-which(is.na(get_dat_ebird$Avibase.ID)==T)
avina<-get_dat_ebird[id,]# unknown 8 sp
avina<-data.frame(spname=avina$spname)

get_dat_ebird_good<-get_dat_ebird%>%filter(!is.na(get_dat_ebird$Avibase.ID))
length(unique(get_dat_ebird_good$spname)) #228 sp found

# now match with ebird,group (if any)
get_dat_ebird_gr<-left_join(avina,avodat_raw,by=c("spname"="eBird.species.group"))
get_dat_ebird_gr_good<-get_dat_ebird_gr%>%filter(!is.na(get_dat_ebird_gr$Avibase.ID))
length(unique(get_dat_ebird_gr_good$spname)) #1 sp found more

id<-which(is.na(get_dat_ebird_gr$Avibase.ID)==T)
avina<-get_dat_ebird_gr[id,]# unknown 7 sp
avina<-data.frame(spname=avina$spname)


# now match with another source
get_dat_BirdTree<-left_join(avina,avodat_raw,by=c("spname"="Species3_BirdTree"))

get_dat_BirdTree_good<-get_dat_BirdTree%>%filter(!is.na(get_dat_BirdTree$Avibase.ID))
length(unique(get_dat_BirdTree_good$spname)) #1 sp found more

id<-which(is.na(get_dat_BirdTree$Avibase.ID)==T)
avina<-get_dat_BirdTree[id,]# unknown 6 sp still
avina<-data.frame(spname=avina$spname)

# now match with another source again
get_dat_BirdLife<-left_join(avina,avodat_raw,by=c("spname"="Species1_BirdLife"))
# no match found!!!!

# so you have to deal with this manually
avina$possible_sp<-NA
avina$comments<-NA
avina2<-data.frame(Avibase.ID=NA)
avina<-cbind(avina2,avina)
write.csv(avina,here("DATA/AVONET/bird_traits_tobe_filledin_min32yr.csv"),row.names=F)

#====================================================================
# first read the manually filled in file
df_manual<-read.csv(here("DATA/AVONET/bird_traits_manually_filled_min32yr.csv"))
df_manual_2<-left_join(df_manual,avodat_raw,by=c("Avibase.ID"="Avibase.ID"))
write.csv(df_manual_2,here("DATA/AVONET/AVONET_data_for_bird_traits_manually_filled_min32yr.csv"),row.names = F)


# ok, now combined all traits
get_dat_ebird_good<-get_dat_ebird_good%>%dplyr::select(-Species1_BirdLife,-eBird.species.group,-Species3_BirdTree)
get_dat_ebird_gr_good<-get_dat_ebird_gr_good%>%dplyr::select(-Species1_BirdLife,-Species2_eBird,-Species3_BirdTree)
get_dat_BirdTree_good<-get_dat_BirdTree_good%>%dplyr::select(-Species1_BirdLife,-Species2_eBird,-eBird.species.group)

birdstraits<-rbind(get_dat_ebird_good,get_dat_ebird_gr_good,get_dat_BirdTree_good)
birdstraits$possible_sp<-birdstraits$spname
birdstraits$comments<-"spname and possible sp are identical, matched sp from AVONET and BBS"
colnames(birdstraits)

df_manual_2<-df_manual_2%>%dplyr::select(colnames(birdstraits))
birdstraits<-rbind(birdstraits,df_manual_2)
birdstraits<-birdstraits[order(birdstraits$spname),]
write.csv(birdstraits,here("DATA/AVONET/bird_traits_from_AVONET_min32yr.csv"),row.names=F)
#birdstraits<-read.csv(here("DATA/AVONET/bird_traits_from_AVONET.csv"))

dftrait<-birdstraits%>%dplyr::select(spname, possible_sp, 
                                 Kipps.Distance, HWI=`Hand-wing.Index`)
dftrait$Kipps.Distance<-as.numeric(dftrait$Kipps.Distance)
dftrait$HWI<-as.numeric(dftrait$HWI)

dftrait2<-dftrait%>%group_by(spname, possible_sp)%>%
  summarize(meanKipps.Distance=mean(Kipps.Distance,na.rm=T),
            meanHWI=mean(HWI,na.rm=T))%>%ungroup()

#=============================== yr_threshold = 32 ==================================================
df32<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_min32yr.csv")) # 32 yrs. = yr_threshold

df_spmeta<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.csv"))
df_spmeta<-df_spmeta%>%dplyr::select(AOU, English_Common_Name, ScientificName)
df32<-left_join(df32,df_spmeta, by="AOU")

mytrait<-left_join(df32, dftrait2, by=c("ScientificName" = "spname"))
write.csv(mytrait,here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_min32yr.csv"),row.names=F)

#=============================== yr_threshold = 36 ==================================================
df36<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_min36yr.csv")) # 36 yrs. = yr_threshold
df_spmeta<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.csv"))
df_spmeta<-df_spmeta%>%dplyr::select(AOU, English_Common_Name, ScientificName)
df36<-left_join(df36,df_spmeta, by="AOU")

mytrait<-left_join(df36, dftrait2, by=c("ScientificName" = "spname"))
write.csv(mytrait,here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_min36yr.csv"),row.names=F)

#=========================
df40<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4.csv")) # 40 yrs. = yr_threshold
df40<-df40%>%dplyr::select(AOU,fLU_ab, fLU_pr,fLU_tas,fLU_pr_avgAprtoAug,fLU_tas_avgAprtoAug, tail)
df_spmeta<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.csv"))
df_spmeta<-df_spmeta%>%dplyr::select(AOU, English_Common_Name, ScientificName)
df40<-left_join(df40,df_spmeta, by="AOU")
mytrait<-left_join(df40, dftrait2, by=c("ScientificName" = "spname"))
write.csv(mytrait,here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits.csv"),row.names=F)

