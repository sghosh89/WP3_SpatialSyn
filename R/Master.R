#-------------------
#path<-dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(path)
#--------------------
library(here)
source(here("R/method_fig.R")) # conceptual figure introducing tail-dependence
source(here("R/get_birds_data_Nyrs_threshold.R")) # wrangle raw data: for 1979-2019 #DONE: 32, 36, 40
source(here("R/prepare_abund_data_Nyrs_threshold.R")) # for a given species, get abundance data in required format #DONE: 32, 36, 40
source(here("R/call_spat_syn_for_abund_Nyrs_threshold.R"))# call for all species' abundance data to compute spat syn #DONE: 32 (running), 36 (running), 40
source(here("R/spat_syn_vs_distance_plot_Nyrs_threshold.R")) # maybe put in suppmat, until 250Km
source(here("R/summary_spat_syn_for_abund_Nyrs_threshold.R")) 

#------------ climate data extraction: remains same ------
source(here("R/download_rawdata_CHELSA.R")) # need a lot of space in your desktop
source(here("DATA/CHELSA_v2/monthly/pr/extract_data_for_uRID_WP3.R"))
source(here("DATA/CHELSA_v2/monthly/tas/extract_data_for_uRID_WP3.R"))

#--------------- spat syn for climate data --------
source(here("R/prepare_climate_data_Nyrs_threshold.R")) # prepare climate data (pr, tas) in required format #DONE: 32, 36, 40
source(here("R/prepare_tas_with_AprtoAug_Nyrs_threshold.R"))# avg data for Apr to Aug #DONE: 32, 36, 40
source(here("R/prepare_pr_with_AprtoAug_Nyrs_threshold.R"))# avg data for Apr to Aug #DONE: 32, 36, 40
source(here("R/call_spat_syn_for_climate_Nyrs_threshold.R")) # compute spat syn for climate data 
source(here("R/summary_spat_syn_for_climate_Nyrs_threshold.R"))# summarize spat_syn for climate data

source(here("R/summarize_res_Nyrs_threshold.R"))# summarize for abundance and climate data

#--------- AVONET: trait data for bird species ----------
# file to get traits for all 78 species 
source(here("R/get_birdtraits_from_AVONET_Nyrs_threshold.R")) 

#-----------------------------------------------------
# Now consider sp. which shows significant tail-dep spatial synchrony

source(here("R/function_to_testsig_abund_Nyrs_threshold.R"))

source(here("R/function_to_testsig_tas_Nyrs_threshold.R"))
source(here("R/function_to_testsig_tas_avgAprtoAug_Nyrs_threshold.R"))

source(here("R/function_to_testsig_pr_Nyrs_threshold.R"))
source(here("R/function_to_testsig_pr_avgAprtoAug_Nyrs_threshold.R"))

source(here("R/distance_sigtaildep_abund_clim_Nyrs_threshold.R"))

source(here("R/visualize_spat_syn_Nyrs_threshold.R"))

#--------- phylogeny data for bird species ----------
source(here("R/get_birdspecies_phylotree_with_Nyrs_threshold.R")) # 180 unique sp. from BirdTREE

source(here("R/phylolm_sig95_0-250km_with_Nyrs_threshold.R"))

source(here("R/phylolm_sig75_0-250km_with_Nyrs_threshold.R"))

source(here("R/phylolm_sig75_0-100km.R"))# for 40yrs
source(here("R/phylolm_sig75_100-250km.R"))# for 40yrs
source(here("R/phylolm_sig95_0-100km.R"))# for 40yrs
source(here("R/phylolm_sig95_100-250km.R"))# for 40yrs

source(here("R/test_phylogenetic_signal_sigtree.R")) # maybe put in suppmat

source(here("R/get_migratory_status.R"))
source(here("R/plot_regression.R"))
source(here("R/plot_td_for_various_months.R"))





