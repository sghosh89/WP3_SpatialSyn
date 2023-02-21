#-------------------
#path<-dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(path)
#--------------------
library(here)
source(here("R/get_birds_data.R")) # wrangle raw data
source(here("R/prepare_abund_data.R")) # for a given species, get abundance data in required format

#source(here("R/compute_spat_syn.R"))# working mode
