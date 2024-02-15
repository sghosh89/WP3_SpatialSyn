# we want to summarize the spat syn across species
library(tidyverse)
library(here)
# first make a table

call_summary_spat_syn_for_climate<-function(df, chosen_rad, climvar,nbin){
  
  summary_spat_syn_for_climate<-c()
  for(i in 1:nrow(df)){
    
    givenAOU<-df$AOU[i]
    
    resloc<-here(paste("RESULTS/AOU_", givenAOU,"/",climvar,"_spatsyn_nbin_",nbin,sep=""))
    
    sdf<-readRDS(paste(resloc,"/summary_df.RDS",sep=""))
    
    inputresloc<-here(paste("RESULTS/AOU_", givenAOU,sep=""))
    distm<-readRDS(here(paste(inputresloc,"/distm_sel.RDS",sep="")))
    
    # summary for a preselected distance category
    
    # now only select sites which are within the interaction radius
    badid_distm_sel<-which(distm<chosen_rad[1] | distm>chosen_rad[2],arr.ind = T)
    distm_good<-distm
    distm_good[badid_distm_sel]<-NA
    
    sdf_selected_range<-sdf
    sdf_selected_range[1,]<-NA
    sdf_selected_range$mindist<-chosen_rad[1]
    sdf_selected_range$maxdist<-chosen_rad[2]
    
    numNA<-colSums(is.na(distm_good))==nrow(distm_good) # count number of columns with all NAs
    sdf_selected_range$nsite<-length(which(numNA==F)) # number of sites with at least one interaction 
    # with other site within chosen radius
    
    good_id<-which(!is.na(distm_good),arr.ind=T)
    sdf_selected_range$nint<-nrow(good_id) # number of pairwise sites within chosen radius
    
    # extract info for selected pairwise sites
    nps<-readRDS(paste(resloc,"/NonParamStat.RDS",sep=""))
    
    #============= tail-dep spatial syn =============
    corlmcoru_good<- nps$Corl - nps$Coru
    corlmcoru_good[badid_distm_sel]<-NA # exclude outside of chosen radius
    
    posnImat_old<-nps$posnI
    posnImat_new<-posnImat_old
    posnImat_new[badid_distm_sel]<-NA
    posnI_new<-which(posnImat_new==1,arr.ind = T)
    
    posnNmat_old<-nps$posnN
    posnNmat_new<-posnNmat_old
    posnNmat_new[badid_distm_sel]<-NA
    posnN_new<-which(posnNmat_new==1,arr.ind = T)
    
    corlmcoru_good[posnI_new]<-NA
    corlmcoru_good[posnN_new]<-NA
    
    #======= now we need to keep only the values that corresponds to corl-coru of abundance =========
    resloc.ab<-here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn_nbin_",nbin,sep=""))
    nps.ab<-readRDS(paste(resloc.ab,"/NonParamStat.RDS",sep=""))
    
    corlmcoru.ab_good<- nps.ab$Corl - nps.ab$Coru
    corlmcoru.ab_good[badid_distm_sel]<-NA # exclude outside of chosen radius
    
    posnImat.ab_old<-nps.ab$posnI
    posnImat.ab_new<-posnImat.ab_old
    posnImat.ab_new[badid_distm_sel]<-NA
    posnI.ab_new<-which(posnImat.ab_new==1,arr.ind = T)
    
    posnNmat.ab_old<-nps.ab$posnN
    posnNmat.ab_new<-posnNmat.ab_old
    posnNmat.ab_new[badid_distm_sel]<-NA
    posnN.ab_new<-which(posnNmat.ab_new==1,arr.ind = T)
    
    corlmcoru.ab_good[posnI.ab_new]<-NA
    corlmcoru.ab_good[posnN.ab_new]<-NA
    
    id.ab.NA<-which(is.na(corlmcoru.ab_good)==T,arr.ind=T)# these are the sites we don't want in climate var too
    corlmcoru_good[id.ab.NA]<-NA
    
    idg<-which(!is.na(corlmcoru_good), arr.ind = T)
    mytab<-data.frame(row=idg[,1],col=idg[,2],CorlmCoru=corlmcoru_good[idg])
    
    if(nrow(mytab)>0){
      saveRDS(mytab, paste(resloc,"/",climvar,"_table_for_sitepair_within_",chosen_rad[1],"_",chosen_rad[2],"_km_nbin_",nbin,".RDS",sep=""))
    }
    #===================================
    
    sdf_selected_range$nL<-sum(corlmcoru_good>0, na.rm = T)
    sdf_selected_range$nU<-sum(corlmcoru_good<0, na.rm = T)
    sdf_selected_range$L<-sum(corlmcoru_good[which(corlmcoru_good>0,arr.ind=T)])
    sdf_selected_range$U<-sum(corlmcoru_good[which(corlmcoru_good<0,arr.ind=T)])
    sdf_selected_range$AOU<-givenAOU
    
    sdf_selected_range<-sdf_selected_range%>%dplyr::select(nsite,nL,nU,L,U,mindist,maxdist,AOU)
    
    saveRDS(sdf_selected_range, paste(resloc,"/summary_df_within_",chosen_rad[1],"_",chosen_rad[2],"_km_nbin_",nbin,".RDS",sep=""))
    
    #print(i)
    #print(sdf_selected_range)
    summary_spat_syn_for_climate<-rbind(summary_spat_syn_for_climate,sdf_selected_range)
  }
  # rearrange
  summary_spat_syn_for_climate<-summary_spat_syn_for_climate%>%select(AOU, mindist, maxdist,
                                                                  nL, nU, L, U)
  
  write.csv(summary_spat_syn_for_climate, here(paste("RESULTS/summary_spat_syn_for_",climvar,"_",
                                                   chosen_rad[1],"_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")),row.names = F)
  
}

# range: 0-250km
df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))
chosen_rad<-c(0,250) # within this distance category

call_summary_spat_syn_for_climate(df=df, chosen_rad=chosen_rad, climvar="pr",nbin=4)
call_summary_spat_syn_for_climate(df=df, chosen_rad=chosen_rad, climvar="pr_avgAprtoAug",nbin=4)


call_summary_spat_syn_for_climate(df=df, chosen_rad=chosen_rad, climvar="tas",nbin=4)

call_summary_spat_syn_for_climate(df=df, chosen_rad=chosen_rad, climvar="tas_avgAprtoAug",nbin=4)
call_summary_spat_syn_for_climate(df=df, chosen_rad=chosen_rad, climvar="tas_avgAprtoJuly",nbin=4)
call_summary_spat_syn_for_climate(df=df, chosen_rad=chosen_rad, climvar="tas_avgMaytoJuly",nbin=4)



