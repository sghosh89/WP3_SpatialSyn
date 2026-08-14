library(here)
library(tidyverse)
library(phylolm)

get_sampling_effort_effect<-function(yr_threshold){
  
  if(yr_threshold==40){
    outputresloc<-here("RESULTS/model_phylolm_sig95_0-250km")
  }else{
    outputresloc<-here(paste0("RESULTS/model_phylolm_sig95_0-250km_min",yr_threshold,"yr"))
  }
  df_sub<-readRDS(here(paste(outputresloc,"/df_sub.RDS",sep="")))
  rownames(df_sub)<-df_sub$Species
  
  ct3<-readRDS(here(paste(outputresloc,"/consensus_tree_with_edgelength.RDS",sep="")))
  
  # Sampling-effort diagnostic among species retained at >=10 sites
  
  df_sub$log_nint <- log(df_sub$nint)
  df_sub$z_log_nint <- as.numeric(scale(df_sub$log_nint))
  
  ct3null<-ct3
  ct3null$node.label<-NULL
  
  set.seed(seed=123)
  # Precipitation
  fitboot1_pr <- phylolm(abs.avg.td.ab.sig~ abs.avg.td.pr5.sig+ HWI + z_log_nint, data=df_sub,
                         phy=ct3,model="lambda", boot=1000)
  saveRDS(fitboot1_pr,here(paste(outputresloc,"model_est_phylolm_boot_with_samplingeffort_P.RDS",sep="/")))
  
  set.seed(seed=123)
  # Temperature
  fitboot1_temp <- phylolm(abs.avg.td.ab.sig~ abs.avg.td.tas5.sig+ HWI + z_log_nint, data=df_sub,
                         phy=ct3,model="lambda", boot=1000)
  saveRDS(fitboot1_temp,here(paste(outputresloc,"model_est_phylolm_boot_with_samplingeffort_T.RDS",sep="/")))
  
  sink(here(paste(outputresloc,"sampling_effort_effect.txt",sep="/")),
       append=TRUE, split=TRUE)
  
  
  cat("============== results using phylolm function with bootstrap: precipitation ============\n")
  print(summary(fitboot1_pr))
  
  cat("============== results using phylolm function with bootstrap: temperature ============\n")
  print(summary(fitboot1_temp))
  
  sink()
  
}

get_sampling_effort_effect(yr_threshold = 32)
get_sampling_effort_effect(yr_threshold = 36)
get_sampling_effort_effect(yr_threshold = 40)
#================================================================

















