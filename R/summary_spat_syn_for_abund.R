# we want to summarize the spat syn across species
library(tidyverse)
library(here)
# first make a table

call_summary_spat_syn_for_abund<-function(df, chosen_rad){
  summary_spat_syn_for_abund<-c()
  for(i in 1:nrow(df)){
    
    givenAOU<-df$AOU[i]
    
    resloc<-here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn",sep=""))
    
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
    
    sdf_selected_range$nind<-nrow(posnI_new)
    sdf_selected_range$nneg<-nrow(posnN_new)
    sdf_selected_range$npos<-sdf_selected_range$nint - (sdf_selected_range$nind + sdf_selected_range$nneg)
    
    corlmcoru_good[posnI_new]<-NA
    corlmcoru_good[posnN_new]<-NA
    
    sdf_selected_range$nL<-sum(corlmcoru_good>0, na.rm = T)
    sdf_selected_range$nU<-sum(corlmcoru_good<0, na.rm = T)
    sdf_selected_range$L<-sum(corlmcoru_good[which(corlmcoru_good>0,arr.ind=T)])
    sdf_selected_range$U<-sum(corlmcoru_good[which(corlmcoru_good<0,arr.ind=T)])
    sdf_selected_range$AOU<-givenAOU
    
    saveRDS(sdf_selected_range, paste(resloc,"/summary_df_within_",chosen_rad[1],"_",chosen_rad[2],"_km.RDS",sep=""))
    
    #print(i)
    #print(sdf_selected_range)
    summary_spat_syn_for_abund<-rbind(summary_spat_syn_for_abund,sdf_selected_range)
  }
  # rearrange
  summary_spat_syn_for_abund<-summary_spat_syn_for_abund%>%select(AOU, mindist, maxdist,
                                                                  nint, nind, npos, nneg, 
                                                                  nL, nU, L, U)
  
  write.csv(summary_spat_syn_for_abund, here(paste("RESULTS/summary_spat_syn_for_abund_",
                                                   chosen_rad[1],"_",chosen_rad[2],"km.csv",sep="")),row.names = F)
  
}

df<-read.csv(here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_w_morethan2sites.csv"))
chosen_rad<-c(0,250) # within this distance category
call_summary_spat_syn_for_abund(df=df, chosen_rad=chosen_rad)

chosen_rad<-c(250,500) # within this distance category
call_summary_spat_syn_for_abund(df=df, chosen_rad=chosen_rad)

chosen_rad<-c(500,750) # within this distance category
call_summary_spat_syn_for_abund(df=df, chosen_rad=chosen_rad)

chosen_rad<-c(750,1000) # within this distance category
call_summary_spat_syn_for_abund(df=df, chosen_rad=chosen_rad)

chosen_rad<-c(0,400) # within this distance category
call_summary_spat_syn_for_abund(df=df, chosen_rad=chosen_rad)

#chosen_rad<-c(1000,3000) # within this distance category
#call_summary_spat_syn_for_abund(df=df, chosen_rad=chosen_rad)

#==================================================



















