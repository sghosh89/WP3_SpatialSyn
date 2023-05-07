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

source(here("R/visualize_spat_syn.R")) # plot synchrony contributions from both tail per diet group

# conceptual figure: AOU_5110
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

# source(here("R/summarize_res.R")) # called in Rmd fifle
# summarizes spatial syn in abund, climate (see source file for details)

#--------- AVONET: traits, bodymass, and phylogeny data for bird species ----------
source(here("R/get_birdtraits_from_AVONET.R"))
source(here("R/get_bodymass_from_AVONET.R"))
# test_hypo code chunk in Rmd file tested hypo related to body traits
source(here("R/get_birdspecies_phylotree.R"))
source(here("R/test_tree.R"))# maybe abundon this later: all random pattern

source(here("R/get_bio_opt_for_species.R")) # optimal P, T

source(here("R/spat_syn_vs_distance_plot.R")) # maybe put as suppmat, also in Rmd

