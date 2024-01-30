# we want to summarize the spat syn across species
library(tidyverse)
library(here)
# first make a table

call_summary_spat_syn_for_abund<-function(df, chosen_rad, nbin){
  summary_spat_syn_for_abund<-c()
  for(i in 1:nrow(df)){
    
    givenAOU<-df$AOU[i]
    
    resloc<-here(paste("RESULTS/AOU_", givenAOU,"/abundance_spatsyn_nbin_",nbin,sep=""))
    
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
    
    sdf_selected_range$nind<-nrow(posnI_new)
    sdf_selected_range$nneg<-nrow(posnN_new)
    sdf_selected_range$npos<-sdf_selected_range$nint - (sdf_selected_range$nind + sdf_selected_range$nneg)
    
    corlmcoru_good[posnI_new]<-NA
    corlmcoru_good[posnN_new]<-NA
    ## keeping only +ve correlation
    #==================================
    
    #========== spatial synchrony: spearman =========
    corval_spear<-nps$corval
    corval_spear[badid_distm_sel]<-NA
    corval_spear[posnI_new]<-NA # keeping all +ve and -ve correlation
    idgs<-which(!is.na(corval_spear), arr.ind = T)
    row.names(idgs)<-NULL
    tempoval<-distm_good[idgs]
    spearval<-corval_spear[idgs]
    
    idgs<-as.data.frame(idgs)
    idgs$dist.KM<-tempoval
    idgs$corval_spear<-spearval
    idgs$row_col<-paste(idgs$row,idgs$col,sep="_")
    #================================================
    #===========================
    idg<-which(!is.na(corlmcoru_good), arr.ind = T)
    mytab<-data.frame(CorlmCoru=corlmcoru_good[idg],
                      row_col=paste(idg[,1],idg[,2],sep="_"))
    idgs<-left_join(idgs,mytab,by="row_col")
    idgs<-idgs%>%dplyr::select(-row_col)
    
    
    if(nrow(idgs)>0){
      saveRDS(idgs, paste(resloc,"/abund_table_for_sitepair_within_",chosen_rad[1],"_",chosen_rad[2],"_km_nbin_",nbin,".RDS",sep=""))
    }
     
    #==========================
    sdf_selected_range$nL<-sum(corlmcoru_good>0, na.rm = T)
    sdf_selected_range$nU<-sum(corlmcoru_good<0, na.rm = T)
    sdf_selected_range$L<-sum(corlmcoru_good[which(corlmcoru_good>0,arr.ind=T)])
    sdf_selected_range$U<-sum(corlmcoru_good[which(corlmcoru_good<0,arr.ind=T)])
    sdf_selected_range$AOU<-givenAOU
    
    saveRDS(sdf_selected_range, paste(resloc,"/summary_df_within_",chosen_rad[1],"_",chosen_rad[2],"_km_nbin_",nbin,".RDS",sep=""))
    
    #print(i)
    #print(sdf_selected_range)
    summary_spat_syn_for_abund<-rbind(summary_spat_syn_for_abund,sdf_selected_range)
  }
  # rearrange
  summary_spat_syn_for_abund<-summary_spat_syn_for_abund%>%select(AOU, mindist, maxdist,
                                                                  nint, nind, npos, nneg, 
                                                                  nL, nU, L, U)
  
  write.csv(summary_spat_syn_for_abund, here(paste("RESULTS/summary_spat_syn_for_abund_",
                                                   chosen_rad[1],"_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")),row.names = F)
  
}

df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))
chosen_rad<-c(0,250) # within this distance category
call_summary_spat_syn_for_abund(df=df, chosen_rad=chosen_rad, nbin=2)
call_summary_spat_syn_for_abund(df=df, chosen_rad=chosen_rad, nbin=4)

#nbin=4
#df1<-read.csv(here(paste("RESULTS/summary_spat_syn_for_abund_",
#                                                 chosen_rad[1],"_",chosen_rad[2],"km_nbin_",nbin,".csv",sep="")))
#df1 has 173 sp. 
#
#id<-which(df1$npos!=0) # 78 sp. were synchronous
#dfsig<-readRDS(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_corlmcoru_sigres_summary.RDS",sep="")))
#dfn0<-dfsig%>%filter(Lsig75!=0 | Usig75!=0)# 59 species showed sig. taildep.
#59/78 = 75% sp. showed tail dep. synchrony (75% Ci scale)
#sum((dfn0$Lsig75 + dfn0$Usig75)>0) # 27 showed LT dep.
#sum((dfn0$Lsig75 + dfn0$Usig75)<0) # 32 showed UT dep.

# 27/78 = 34% sp. showed LT dep.
# 32/78 = 41% sp. showed UT dep.



