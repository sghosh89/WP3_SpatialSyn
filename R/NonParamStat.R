# ------------------------------------------------------------------------------------------------
# THIS CODE IS WRITTEN TO TEST NON-PARAMETRIC STATS based on partial Spearmn correlation (Cor stat) 
#              FOR ANY TWO POSITIVELY CORRELATED SPECIES PAIR
#--------------------------------------------------------------------------------------------------
library(here)
source(here("R/CorlCoru.R"))
source(here("R/vivj_matrix.R"))
#---------------------------------------------------------------
#Processing function for a copula approach to synchrony
#
#---Input---------
# m : output matrix from vivj_matrix.R
# nbin : number of bins used to get npa stats

#------Output - A list with these elements-------------------------------
# ranks :       A dataframe with 2 columns, one for 
#               ds1 and one for ds2, corresponding to samples from the copula.
# spear   :     Spearman correlation (single number)
# kend      :   Kendall correlation (single number)
# Corl,Coru  :  covariance based statistics

copsync<-function(m,nbin){
  
  vi<-m[,1]
  vj<-m[,2]
  
  #get mean and variance
  vi_mean<-mean(vi)
  vj_mean<-mean(vj)
  var_vi<-var(vi)
  var_vj<-var(vj)
  
  spear<-NA
  kend<-NA
  Corl<-NA
  Coru<-NA
  
  
  if (length(vi)>0){
    #get spear
    spear<-cor(vi,vj, method ="spearman") 
    
    #get kend
    kend<-cor(vi, vj, method ="kendall") 
    
    #----------------------------------------------------STATISTICS :2 ---------------------------------------------------------
    #get Corl, Coru   (covariance based new stat)
    
    # if(npa_stats=="cor"){
    stat2<-CorlCoru(vi,vj,nbin=nbin)
    Corl<-stat2[1]
    Coru<-stat2[2]
    # }
  }
  return(list(ranks=data.frame(Rki=vi,Rkj=vj),
              spear=spear,kend=kend,
              Corl=Corl,Coru=Coru))
}

#---------------------------------------------------------------------------------------------------------------------------

#Calling the above on all pairs of several time series and returning npa stats
#-------------Input---------------------------------
# d_allsite  :    A list of data frames, each with columns Year and Dat
#               The years are assumed to be sequential and all included,
#               though there may be NAs in Dat and the years may not
#               be all the same for ds1 and ds2.
# pfname   :    Filename (without extension) prepended to plot files saved.
# nbin : number of bins used to compute npa stats

#----------Output - A list of length 1 with these elements-----------------

# spear            A matrix of spearman results, length(d) by length(d) (significantly -ve correlated cells 
#                    had +ve correlation, 
#                     as we take one of timeseries reversed but can track with posnN indices)

# kend             A matrix of kendall results, length(d) by length(d)

# Corl             A matrix of Cl results, length(d) by length(d)
# Coru             A matrix of Cu results, length(d) by length(d)

#posnI             A matrix for indices of indep. sp. pair
#posnN             A matrix for indices of negatively correlated sp. pair

#corval           A matrix of spearman results but differ with spear in the sense that it contains all +ve or -ve Spearman
#                 correlation values, whether significant or not. And also the true correlation. (same as 'corval' output 
#                                 from vivj_matrix.R)

multcall<-function(d_allsite,resloc,nbin){
  
  lensite<-length(d_allsite)
  
  s <- strsplit(names(d_allsite), " ");
  cnms<-names(d_allsite)
  
  #first initialize result receptacles for the output
  spear<-matrix(NA,lensite,lensite)
  colnames(spear) <- cnms
  rownames(spear) <-names(d_allsite)
  
  kend<-matrix(NA,lensite,lensite)
  colnames(kend) <- colnames(spear)
  rownames(kend) <-rownames(spear)
  
  Corl<-matrix(NA,lensite,lensite)
  colnames(Corl) <- colnames(spear)
  rownames(Corl) <-rownames(spear)
  
  Coru<-matrix(NA,lensite,lensite)
  colnames(Coru) <- colnames(spear)
  rownames(Coru) <-rownames(spear)
  
  corval<-spear
  pval_BiCopIndep<-spear
  posnI<-spear
  posnN<-spear
  
  #------------------- PLOT :  copula_for_all_sp pair ----------------
  pdf(paste(resloc,"AllCops_nbin_",nbin,".pdf",sep=""),width=6*lensite, height=6*lensite)
  op<-par(mfrow=c(lensite,lensite),mar=c(3,3,3,3), mgp=c(1.5,0.5,0))
  
  for (i in c(1:lensite)){
    for (j in c(1:lensite)){
      
      #if(i!=j){
      ms<-vivj_matrix(d_allsite=d_allsite,i=i,j=j,level=0.05,ploton=T)
      m<-ms$mat
      
      corval[i,j]<-ms$corval
      pval_BiCopIndep[i,j]<-ms$IndepTestRes
      
      thisres<-copsync(m=m,nbin=nbin)
      
      spear[i,j]<-thisres$spear
      kend[i,j]<-thisres$kend
      
      Corl[i,j]<-thisres$Corl
      Coru[i,j]<-thisres$Coru
      
      #plot(thisres$ranks$Rki,thisres$ranks$Rkj,type='p',col=rgb(0,0,0,.2),pch=19,xlim=c(0,1),
      #     ylim=c(0,1),xlab=expression(u[i]),ylab=expression(v[j]),cex.lab=2)
      #mtext(paste0("[ i, j ] ="," [",i,",",j,"] ", ","," n=",dim(thisres$ranks)[1]),
      #      side = 3, line=0.15, adj=0.5, col="red")
      #}
    }
  }
  par(op)
  dev.off()
  
  #-------------------------------------------------------------------------
  
  
   level<-0.05
   posnIind<-which(pval_BiCopIndep>=level, arr.ind = T) #indices of indep. pair
   posnNind<-which(pval_BiCopIndep<level & corval <0, arr.ind = T) #indices of significantly neg. correlated pair
 
  
  posnI[posnIind]<-1
  posnN[posnNind]<-1
  
  #-------- saving only lower triangular part of the result matrix ------
  spear[upper.tri(spear,diag=T)]<-NA
  kend[upper.tri(kend,diag=T)]<-NA
  Corl[upper.tri(Corl,diag=T)]<-NA
  Coru[upper.tri(Coru,diag=T)]<-NA
  corval[upper.tri(corval,diag=T)]<-NA
  posnI[upper.tri(posnI,diag=T)]<-NA
  posnN[upper.tri(posnN,diag=T)]<-NA
  #---------------------------------------------
  
  #if(npa_stats=="cor"){
  res<-list(spear=spear,kend=kend,
            Corl=Corl,Coru=Coru,
            posnI=posnI,
            posnN=posnN,
            corval=corval)
  #}
  
  return(res)
} 






