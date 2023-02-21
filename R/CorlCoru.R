#-----------------------------------------------------------------------------------------------------------------------
# THIS CODE CONTAINS STATISTICS FOR LOOKING AT TAIL DEPENDENCE
#-----------------------------------------------------------------------------------------------------------------------
library(here)
source(here("R/Corbds.R"))

#--------------------------- correlation based Stat ---------------------------------------------------------

CorlCoru<-function(vi,vj,nbin){
  
  lb<-0
  ub<-1/nbin
  
  Corl<-Corbds(vi,vj,lb=lb,ub=ub)
  Coru<-Corbds(vi,vj,lb=1-ub,ub=1-lb)
  
  return(c(Corl,Coru))
}