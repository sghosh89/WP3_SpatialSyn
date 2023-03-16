#-------------------
#path<-dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(path)
#--------------------
library(here)
source(here("R/get_birds_data.R")) # wrangle raw data
source(here("R/prepare_abund_data.R")) # for a given species, get abundance data in required format
source(here("R/call_spat_syn_for_abund.R"))# call for all species' abundance data to compute spat syn
source(here("R/summary_spat_syn_for_abund.R")) # summarize birds' abundance spat syn result
source(here("R/diet_cat.R")) # get diet and foraging info for each bird species

# even with IUCN status, need to edit below plot a bit 
#source(here("R/plot_summary_dietcat.R")) # plot synchrony contributions from both tail per diet group

source(here("R/visualize_spat_syn.R")) # plot synchrony contributions from both tail per diet group

#------------ climate data extraction ------
source(here("R/download_rawdata_CHELSA.R")) # need a lot of space in your desktop
source(here("DATA/CHELSA_v2/monthly/pr/extract_data_for_uRID_WP3.R"))
source(here("DATA/CHELSA_v2/monthly/tas/extract_data_for_uRID_WP3.R"))# fix this
source(here("DATA/CHELSA_v2/monthly/tasmax/extract_data_for_uRID_WP3.R"))
source(here("DATA/CHELSA_v2/monthly/tasmin/extract_data_for_uRID_WP3.R"))

#--------------- spat syn for climate data --------
source(here("R/prepare_climate_data.R")) # prepare climate data (pr, tas, tasmax, tasmin) in required format
source(here("R/call_spat_syn_for_climate.R")) # compute spat syn for climate data (pr, tas, tasmax, tasmin) 
source(here("R/summary_spat_syn_for_climate.R"))# summarize spat_syn for climate data






