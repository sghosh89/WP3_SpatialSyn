library(here)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
# Now, we want to know the migratory status for the 78 bird species
df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin_min32yr.csv"))
df<-df%>%dplyr::select(AOU,English_Common_Name,ScientificName,BirdTreeName)

migdf<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/MigrantNonBreeder/MigrantSummary.csv"))
migNB_AOU<-unique(migdf$AOU)

id<-which(df$AOU %in% migNB_AOU) # 112 found migrant nonbreeder
idNA<-which(df$AOU %notin% migNB_AOU) #24 sp. not found

df$mig_NB_BBS<-NA
df$mig_NB_BBS[id]<-1

df$comments<-NA
df$source<-NA
write.csv(df,here("RESULTS/species_0_250km_nbin_4_tailsig95_migstatus_tobefilled_min32yr.csv"),row.names = F)
#====================================
dfmig<-read.csv(here("RESULTS/species_0_250km_nbin_4_tailsig95_migstatus_manuallyfilled_min32yr.csv"))
table(df$mig_NB_BBS)
dfmig<-dfmig%>%dplyr::select(AOU,mig_NB_BBS)
#===========================
df<-read.csv(here("DATA/BirdTree/species_0_250km_nbin_4_filledin_min32yr.csv"))
df$newBT<-gsub(" ", "_", df$BirdTreeName)
df<-df%>%dplyr::select(AOU,newBT,ScientificName,BirdTreeName)

nbin<-4
dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km_min32yr.csv",sep="")))
dfsig<-left_join(dfsig,df,by="AOU")


#dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits.csv"))
dft<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4_with_speciestraits_min32yr.csv"))
dft<-dft%>%dplyr::select(ScientificName,kipps=meanKipps.Distance,HWI=meanHWI)
dfsig<-left_join(dfsig,dft,by="ScientificName")

dfsig_mig<-left_join(dfsig,dfmig,by="AOU")

dfsig_migUT<-dfsig_mig%>%filter(tail95=="UT")
as.data.frame(table(dfsig_migUT$mig_NB_BBS)) # 82 out of 88 sp. are migratory

dfsig_migLT<-dfsig_mig%>%filter(tail95=="LT")
as.data.frame(table(dfsig_migLT$mig_NB_BBS)) # 39 out of 44 sp. are migratory

as.data.frame(table(dfsig_mig$mig_NB_BBS))
# in total 121 out of 132 species are migratory






