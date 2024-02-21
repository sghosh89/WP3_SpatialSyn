# write a function to get summarized value of sig tail dep (75%CI) within a given distance boundary
library(here)
library(tidyverse)

distance_sigtaildep_abund_clim<-function(df,nbin=4,target_dist_cat){
  
  #nbin<-4
  #df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4.csv"))
  df<-df%>%dplyr::select(AOU)
  
  df$Lsig75ab<-NA
  df$Usig75ab<-NA
  df$Lsig95ab<-NA
  df$Usig95ab<-NA
  
  df$Lsig75tas<-NA
  df$Usig75tas<-NA
  df$Lsig95tas<-NA
  df$Usig95tas<-NA
  
  df$Lsig75tas5<-NA
  df$Usig75tas5<-NA
  df$Lsig95tas5<-NA
  df$Usig95tas5<-NA
  
  df$Lsig75pr<-NA
  df$Usig75pr<-NA
  df$Lsig95pr<-NA
  df$Usig95pr<-NA
  
  df$Lsig75pr5<-NA
  df$Usig75pr5<-NA
  df$Lsig95pr5<-NA
  df$Usig95pr5<-NA
  
  
  #target_dist_cat<-c(0,100)
  
  for(i in 1:nrow(df)){
    
    givenAOU<-df$AOU[i]
    
    # abundance
    tempoab<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
    #tas
    tempotas<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/tas_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
    #tas_avgAprtoAug
    tempotas5<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/tas_avgAprtoAug_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
    #tas_avgAprtoJuly
    #tempotas4<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/tas_avgAprtoJuly_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
    #tas_avgMaytoJuly
    #tempotas3<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/tas_avgMaytoJuly_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
    #pr
    tempopr<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/pr_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
    #pr_avgAprtoAug
    tempopr5<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/pr_avgAprtoAug_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
     
    #============= first for abundance ==================
    tempoab<-tempoab%>%filter(dist.KM>target_dist_cat[1] & dist.KM<=target_dist_cat[2])
    
    tempo75<-tempoab%>%filter(sig75==1)
    
    if(nrow(tempo75)==0){
      Lsig75<-Usig75<-0
    }else{
      Lsig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual>0)])
      Usig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual<0)])
    }
    
    tempo95<-tempoab%>%filter(sig95==1)
    
    if(nrow(tempo95)==0){
      Lsig95<-Usig95<-0
    }else{
      Lsig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual>0)])
      Usig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual<0)])
    }
    
    df$Lsig75ab[i]<-Lsig75
    df$Usig75ab[i]<-Usig75
    df$Lsig95ab[i]<-Lsig95
    df$Usig95ab[i]<-Usig95
    
    #============= now for tas ==================
    
    tempotas<-tempotas%>%filter(dist.KM>target_dist_cat[1] & dist.KM<=target_dist_cat[2])
    
    tempo75<-tempotas%>%filter(sig75==1)
    
    if(nrow(tempo75)==0){
      Lsig75<-Usig75<-0
    }else{
      Lsig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual>0)])
      Usig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual<0)])
    }
    
    tempo95<-tempotas%>%filter(sig95==1)
    
    if(nrow(tempo95)==0){
      Lsig95<-Usig95<-0
    }else{
      Lsig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual>0)])
      Usig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual<0)])
    }
    
    df$Lsig75tas[i]<-Lsig75
    df$Usig75tas[i]<-Usig75
    df$Lsig95tas[i]<-Lsig95
    df$Usig95tas[i]<-Usig95
    
    #============= now for tas5 ==================
    
    tempotas5<-tempotas5%>%filter(dist.KM>target_dist_cat[1] & dist.KM<=target_dist_cat[2])
    
    tempo75<-tempotas5%>%filter(sig75==1)
    
    if(nrow(tempo75)==0){
      Lsig75<-Usig75<-0
    }else{
      Lsig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual>0)])
      Usig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual<0)])
    }
    
    tempo95<-tempotas5%>%filter(sig95==1)
    
    if(nrow(tempo95)==0){
      Lsig95<-Usig95<-0
    }else{
      Lsig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual>0)])
      Usig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual<0)])
    }
    
    df$Lsig75tas5[i]<-Lsig75
    df$Usig75tas5[i]<-Usig75
    df$Lsig95tas5[i]<-Lsig95
    df$Usig95tas5[i]<-Usig95
    
    #============= now for pr ==================
    
    tempopr<-tempopr%>%filter(dist.KM>target_dist_cat[1] & dist.KM<=target_dist_cat[2])
    
    tempo75<-tempopr%>%filter(sig75==1)
    
    if(nrow(tempo75)==0){
      Lsig75<-Usig75<-0
    }else{
      Lsig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual>0)])
      Usig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual<0)])
    }
    
    tempo95<-tempopr%>%filter(sig95==1)
    
    if(nrow(tempo95)==0){
      Lsig95<-Usig95<-0
    }else{
      Lsig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual>0)])
      Usig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual<0)])
    }
    
    df$Lsig75pr[i]<-Lsig75
    df$Usig75pr[i]<-Usig75
    df$Lsig95pr[i]<-Lsig95
    df$Usig95pr[i]<-Usig95
    
    #============= now for pr: avg over 5 months: Apr to Aug ==================
    
    tempopr5<-tempopr5%>%filter(dist.KM>target_dist_cat[1] & dist.KM<=target_dist_cat[2])
    
    tempo75<-tempopr5%>%filter(sig75==1)
    
    if(nrow(tempo75)==0){
      Lsig75<-Usig75<-0
    }else{
      Lsig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual>0)])
      Usig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual<0)])
    }
    
    tempo95<-tempopr5%>%filter(sig95==1)
    
    if(nrow(tempo95)==0){
      Lsig95<-Usig95<-0
    }else{
      Lsig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual>0)])
      Usig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual<0)])
    }
    
    df$Lsig75pr5[i]<-Lsig75
    df$Usig75pr5[i]<-Usig75
    df$Lsig95pr5[i]<-Lsig95
    df$Usig95pr5[i]<-Usig95
    
    print(i)
  }
  
  saveRDS(df,here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,
                        "_corlmcoru_sigres_summary_",target_dist_cat[1],"-",
                        target_dist_cat[2],"Km.RDS",sep="")))
  
}

nbin<-4
df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4.csv"))
target_dist_cat<-c(0,100)
distance_sigtaildep_abund_clim(df=df,nbin=nbin,target_dist_cat=target_dist_cat)

nbin<-4
df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4.csv"))
target_dist_cat<-c(100,250)
distance_sigtaildep_abund_clim(df=df,nbin=nbin,target_dist_cat=target_dist_cat)

nbin<-4
df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4.csv"))
target_dist_cat<-c(0,250)
distance_sigtaildep_abund_clim(df=df,nbin=nbin,target_dist_cat=target_dist_cat)

#===================================

#===================================

get_summary_csv<-function(nbin=4,target_dist_cat,siglevel=75){
  
  dff<-readRDS(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,
                          "_corlmcoru_sigres_summary_",target_dist_cat[1],"-",
                          target_dist_cat[2],"Km.RDS",sep="")))
  
  if(siglevel==75){
    dff<-dff%>%filter(Lsig75ab!=0 | Usig75ab!=0)# 59 species having sig taildep in abundance
    
    dff$fab.sig<-(dff$Lsig75ab+dff$Usig75ab)/(dff$Lsig75ab+abs(dff$Usig75ab))
    dff$ftas.sig<-(dff$Lsig75tas+dff$Usig75tas)/(dff$Lsig75tas+abs(dff$Usig75tas))
    
    dff$ftas5.sig<-(dff$Lsig75tas5+dff$Usig75tas5)/(dff$Lsig75tas5+abs(dff$Usig75tas5))
    #dff$ftas4.sig<-(dff$Lsig75tas4+dff$Usig75tas4)/(dff$Lsig75tas4+abs(dff$Usig75tas4))
    #dff$ftas3.sig<-(dff$Lsig75tas3+dff$Usig75tas3)/(dff$Lsig75tas3+abs(dff$Usig75tas3))
    
    dff$fpr.sig<-(dff$Lsig75pr+dff$Usig75pr)/(dff$Lsig75pr+abs(dff$Usig75pr))
    dff$fpr5.sig<-(dff$Lsig75pr5+dff$Usig75pr5)/(dff$Lsig75pr5+abs(dff$Usig75pr5))
    
    
    tail75<-(dff$Lsig75ab+dff$Usig75ab)
    dff$tail75<-ifelse(tail75<0,"UT","LT")
    dff$tail75<-as.factor(dff$tail75)
    
    dff$abs.tot.td.ab.sig<-dff$Lsig75ab+abs(dff$Usig75ab)
    dff$abs.tot.td.tas5.sig<-dff$Lsig75tas5+abs(dff$Usig75tas5)
    dff$abs.tot.td.pr5.sig<-dff$Lsig75pr5+abs(dff$Usig75pr5)
    
    dff$tot.td.ab.sig<-dff$Lsig75ab+dff$Usig75ab
    dff$tot.td.tas5.sig<-dff$Lsig75tas5+dff$Usig75tas5
    dff$tot.td.pr5.sig<-dff$Lsig75pr5+dff$Usig75pr5
    
    dff<-dff%>%dplyr::select(AOU,tail75,
                             abs.tot.td.ab.sig,
                             tot.td.ab.sig,
                             fab.sig,
                             ftas.sig,
                             ftas5.sig,
                             abs.tot.td.tas5.sig,
                             tot.td.tas5.sig,
                             #ftas4.sig,
                             #ftas3.sig,
                             fpr.sig,
                             fpr5.sig,
                             abs.tot.td.pr5.sig,
                             tot.td.pr5.sig)
    write.csv(dff,here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail75sig_summary_",
                             target_dist_cat[1],"-",
                             target_dist_cat[2],"Km.csv",sep="")), row.names = F)
  }
  
  if(siglevel==95){
    dff<-dff%>%filter(Lsig95ab!=0 | Usig95ab!=0)
    
    dff$fab.sig<-(dff$Lsig95ab+dff$Usig95ab)/(dff$Lsig95ab+abs(dff$Usig95ab))
    dff$ftas.sig<-(dff$Lsig95tas+dff$Usig95tas)/(dff$Lsig95tas+abs(dff$Usig95tas))
    
    dff$ftas5.sig<-(dff$Lsig95tas5+dff$Usig95tas5)/(dff$Lsig95tas5+abs(dff$Usig95tas5))
    #dff$ftas4.sig<-(dff$Lsig95tas4+dff$Usig95tas4)/(dff$Lsig95tas4+abs(dff$Usig95tas4))
    #dff$ftas3.sig<-(dff$Lsig95tas3+dff$Usig95tas3)/(dff$Lsig95tas3+abs(dff$Usig95tas3))
    
    dff$fpr.sig<-(dff$Lsig95pr+dff$Usig95pr)/(dff$Lsig95pr+abs(dff$Usig95pr))
    dff$fpr5.sig<-(dff$Lsig95pr5+dff$Usig95pr5)/(dff$Lsig95pr5+abs(dff$Usig95pr5))
    
    tail95<-(dff$Lsig95ab+dff$Usig95ab)
    dff$tail95<-ifelse(tail95<0,"UT","LT")
    dff$tail95<-as.factor(dff$tail95)
    
    dff$abs.tot.td.ab.sig<-dff$Lsig95ab+abs(dff$Usig95ab)
    dff$abs.tot.td.tas5.sig<-dff$Lsig95tas5+abs(dff$Usig95tas5)
    dff$abs.tot.td.pr5.sig<-dff$Lsig95pr5+abs(dff$Usig95pr5)
    
    dff$tot.td.ab.sig<-dff$Lsig95ab+dff$Usig95ab
    dff$tot.td.tas5.sig<-dff$Lsig95tas5+dff$Usig95tas5
    dff$tot.td.pr5.sig<-dff$Lsig95pr5+dff$Usig95pr5
    
    dff<-dff%>%dplyr::select(AOU,tail95,
                             abs.tot.td.ab.sig,
                             tot.td.ab.sig,
                             fab.sig,
                             ftas.sig,
                             ftas5.sig,
                             abs.tot.td.tas5.sig,
                             tot.td.tas5.sig,
                             #ftas4.sig,
                             #ftas3.sig,
                             fpr.sig,
                             fpr5.sig,
                             abs.tot.td.pr5.sig,
                             tot.td.pr5.sig)
    
    write.csv(dff,here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_",
                             target_dist_cat[1],"-",
                             target_dist_cat[2],"Km.csv",sep="")), row.names = F)
  }
  
} 

# call the above function
nbin<-4
target_dist_cat<-c(0,100)
siglevel=75
get_summary_csv(nbin=nbin,target_dist_cat = target_dist_cat, siglevel=siglevel)

nbin<-4
target_dist_cat<-c(100,250)
siglevel=75
get_summary_csv(nbin=nbin,target_dist_cat = target_dist_cat, siglevel=siglevel)

nbin<-4
target_dist_cat<-c(0,250)
siglevel=75
get_summary_csv(nbin=nbin,target_dist_cat = target_dist_cat, siglevel=siglevel)



