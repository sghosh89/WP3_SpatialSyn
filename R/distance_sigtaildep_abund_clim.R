# write a function to get summarized value of sig tail dep (75%CI) within a given distance boundary


distance_sigtaildep_abund_clim<-function(df,nbin=4,target_dist_cat){
  
  #nbin<-4
  #df<-read.csv(here("RESULTS/df_abund_climate_spatsyn_0_250km_nbin_4.csv"))
  df<-df%>%dplyr::select(AOU)
  
  df$Lsig75ab<-NA
  df$Usig75ab<-NA
  df$Lsig95ab<-NA
  df$Usig95ab<-NA
  
  df$Lsig75tasmax<-NA
  df$Usig75tasmax<-NA
  df$Lsig95tasmax<-NA
  df$Usig95tasmax<-NA
  
  df$Lsig75tasmax5<-NA
  df$Usig75tasmax5<-NA
  df$Lsig95tasmax5<-NA
  df$Usig95tasmax5<-NA
  
  #target_dist_cat<-c(0,100)
  
  for(i in 1:nrow(df)){
    
    givenAOU<-df$AOU[i]
    
    # abundance
    tempoab<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
    #tasmax
    tempotasmax<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/tasmax_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
    #tasmax_avgAprtoAug
    tempotasmax5<-readRDS(here(paste("RESULTS/AOU_", givenAOU,"/tasmax_avgAprtoAug_spatsyn_nbin_",nbin,"/corlmcoru_sigres.RDS",sep="")))
    
    # check
    #ch1<-any(tempoab$dist.KM==tempotasmax$dist.KM)==F
    #ch2<-any(tempoab$dist.KM==tempotasmax5$dist.KM)==F
    #ch3<-any(ch1,ch2)==T
    
    #if(ch3==T){
    # print("something went wrong!!!")
    #}
    
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
    
    #============= now for tasmax ==================
    
    tempotasmax<-tempotasmax%>%filter(dist.KM>target_dist_cat[1] & dist.KM<=target_dist_cat[2])
    
    tempo75<-tempotasmax%>%filter(sig75==1)
    
    if(nrow(tempo75)==0){
      Lsig75<-Usig75<-0
    }else{
      Lsig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual>0)])
      Usig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual<0)])
    }
    
    tempo95<-tempotasmax%>%filter(sig95==1)
    
    if(nrow(tempo95)==0){
      Lsig95<-Usig95<-0
    }else{
      Lsig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual>0)])
      Usig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual<0)])
    }
    
    df$Lsig75tasmax[i]<-Lsig75
    df$Usig75tasmax[i]<-Usig75
    df$Lsig95tasmax[i]<-Lsig95
    df$Usig95tasmax[i]<-Usig95
    
    #============= now for tasmax: avg over 5 months: Apr to Aug ==================
    
    tempotasmax5<-tempotasmax5%>%filter(dist.KM>target_dist_cat[1] & dist.KM<=target_dist_cat[2])
    
    tempo75<-tempotasmax5%>%filter(sig75==1)
    
    if(nrow(tempo75)==0){
      Lsig75<-Usig75<-0
    }else{
      Lsig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual>0)])
      Usig75<-sum(tempo75$corlmcoru_actual[which(tempo75$corlmcoru_actual<0)])
    }
    
    tempo95<-tempotasmax5%>%filter(sig95==1)
    
    if(nrow(tempo95)==0){
      Lsig95<-Usig95<-0
    }else{
      Lsig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual>0)])
      Usig95<-sum(tempo95$corlmcoru_actual[which(tempo95$corlmcoru_actual<0)])
    }
    
    df$Lsig75tasmax5[i]<-Lsig75
    df$Usig75tasmax5[i]<-Usig75
    df$Lsig95tasmax5[i]<-Lsig95
    df$Usig95tasmax5[i]<-Usig95
    
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
    dff$ftasmax.sig<-(dff$Lsig75tasmax+dff$Usig75tasmax)/(dff$Lsig75tasmax+abs(dff$Usig75tasmax))
    dff$ftasmax5.sig<-(dff$Lsig75tasmax5+dff$Usig75tasmax5)/(dff$Lsig75tasmax5+abs(dff$Usig75tasmax5))
    
    
    tail75<-(dff$Lsig75ab+dff$Usig75ab)
    dff$tail75<-ifelse(tail75<0,"UT","LT")
    dff$tail75<-as.factor(dff$tail75)
    
    dff<-dff%>%dplyr::select(AOU,tail75,fab.sig,ftasmax.sig,ftasmax5.sig)
    write.csv(dff,here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail75sig_summary_",
                             target_dist_cat[1],"-",
                             target_dist_cat[2],"Km.csv",sep="")), row.names = F)
  }
  
  if(siglevel==95){
    dff<-dff%>%filter(Lsig95ab!=0 | Usig95ab!=0)
    
    dff$fab.sig<-(dff$Lsig95ab+dff$Usig95ab)/(dff$Lsig95ab+abs(dff$Usig95ab))
    dff$ftasmax.sig<-(dff$Lsig95tasmax+dff$Usig95tasmax)/(dff$Lsig95tasmax+abs(dff$Usig95tasmax))
    dff$ftasmax5.sig<-(dff$Lsig95tasmax5+dff$Usig95tasmax5)/(dff$Lsig95tasmax5+abs(dff$Usig95tasmax5))
    
    
    tail95<-(dff$Lsig95ab+dff$Usig95ab)
    dff$tail95<-ifelse(tail95<0,"UT","LT")
    dff$tail95<-as.factor(dff$tail95)
    
    dff<-dff%>%dplyr::select(AOU,tail95,fab.sig,ftasmax.sig,ftasmax5.sig)
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
