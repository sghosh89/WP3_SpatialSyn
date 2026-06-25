
#rm(list=ls())
library(here)
library(dplyr)


yr_threshold <- 32 
target_dist_cat <- c(0, 250) 
nbin <- 4 
siglevel <- "95%CI"

if(yr_threshold==40){
  df1<-readRDS(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,
                          "_corlmcoru_sigres_summary_",target_dist_cat[1],"-",
                          target_dist_cat[2],"Km.RDS",sep="")))
}else{
  df1<-readRDS(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,
                          "_corlmcoru_sigres_summary_",target_dist_cat[1],"-",
                          target_dist_cat[2],"Km_min",yr_threshold,"yr.RDS",sep="")))
}

if(siglevel=="95%CI"){
  df1<-df1%>%filter(Lsig95ab!=0 | Usig95ab!=0)
  df1<-df1%>%dplyr::select(AOU,nint,Lsigab=Lsig95ab, Usigab=Usig95ab)
}

df_spmeta<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.csv"))
df_spmeta<-df_spmeta%>%dplyr::select(AOU, English_Common_Name, ScientificName)
df<-left_join(df1,df_spmeta, by="AOU")

df<-df%>%select(AOU, nint, English_Common_Name, ScientificName)%>%filter(nint>=45)%>%select(-nint) # show this table into suppmat
write.csv(df, here("Results/SpeciesList_0_250km_nbin_4_min32yr_min10sites.csv"),row.names=F) # 124 species
