library(here)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
# Now, we want to know the migratory status for the 78 bird species
df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin.csv"))
df<-df%>%dplyr::select(AOU,English_Common_Name,ScientificName,BirdTreeName)

migdf<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/MigrantNonBreeder/MigrantSummary.csv"))
migNB_AOU<-unique(migdf$AOU)

id<-which(df$AOU %in% migNB_AOU) # 62 found migrant nonbreeder
idNA<-which(df$AOU %notin% migNB_AOU) #16 sp. not found

df$mig_NB_BBS<-NA
df$mig_NB_BBS[id]<-1

df$comments<-NA
df$source<-NA
write.csv(df,here("RESULTS/species_0_250km_nbin_4_tailsig75_migstatus_tobefilled.csv"),row.names = F)
#====================================
dfmig<-read.csv(here("RESULTS/species_0_250km_nbin_4_tailsig75_migstatus_manuallyfilled.csv"))
table(df$mig_NB_BBS)
dfmig<-dfmig%>%dplyr::select(AOU,mig_NB_BBS)
#===========================
df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin_min32yr.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
df<-df%>%dplyr::select(AOU,newBT,ScientificName,BirdTreeName)

nbin<-4
dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail75sig_summary_0-250Km.csv",sep="")))
dfsig<-left_join(dfsig,df,by="AOU")


dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits.csv"))
dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
dfsig<-left_join(dfsig,dft,by="ScientificName")

dfsig_mig<-left_join(dfsig,dfmig,by="AOU")

dfsig_migUT<-dfsig_mig%>%filter(tail75=="UT")
table(dfsig_migUT$mig_NB_BBS) # 30 out of 32 sp. are migratory

dfsig_migLT<-dfsig_mig%>%filter(tail75=="LT")
table(dfsig_migLT$mig_NB_BBS) # 22 out of 27 sp. are migratory

# in total 52 out of 59 species are migratory






