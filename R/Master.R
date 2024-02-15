#-------------------
#path<-dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(path)
#--------------------
library(here)
source(here("R/method_fig.R")) # conceptual figure introducing tail-dependence
source(here("R/get_birds_data.R")) # wrangle raw data: for 1979-2019
source(here("R/prepare_abund_data.R")) # for a given species, get abundance data in required format
source(here("R/call_spat_syn_for_abund.R"))# call for all species' abundance data to compute spat syn
source(here("R/spat_syn_vs_distance_plot.R")) # maybe put in suppmat, until 250Km
source(here("R/summary_spat_syn_for_abund.R")) # summarize birds' abundance spat syn result (373 sp within 0-250 km distance with more than 2 sites)
source(here("R/diet_cat.R")) # get diet and foraging info for each bird species

#source(here("R/visualize_spat_syn.R")) # plot synchrony contributions from both tail per diet group between 0-250 Km between sites distance, total 263 sp. there 
# later visualize with significant tail dep result only

# conceptual figure: AOU_5110
#------------ climate data extraction ------
source(here("R/download_rawdata_CHELSA.R")) # need a lot of space in your desktop
source(here("DATA/CHELSA_v2/monthly/pr/extract_data_for_uRID_WP3.R"))
source(here("DATA/CHELSA_v2/monthly/tas/extract_data_for_uRID_WP3.R"))
source(here("DATA/CHELSA_v2/monthly/tasmax/extract_data_for_uRID_WP3.R"))

#--------------- spat syn for climate data --------
source(here("R/prepare_climate_data.R")) # prepare climate data (pr, tas) in required format

source(here("R/prepare_tas_with_AprtoAug.R"))# avg data for Apr to Aug
source(here("R/prepare_pr_with_AprtoAug.R"))# avg data for Apr to Aug

source(here("R/prepare_tas_with_AprtoJuly.R"))# avg data for Apr to July
source(here("R/prepare_tas_with_MaytoJuly.R"))# avg data for May to July

source(here("R/call_spat_syn_for_climate.R")) # compute spat syn for climate data 

source(here("R/summary_spat_syn_for_climate.R"))# summarize spat_syn for climate data

# ftd_abund vs ftd_climate plot for both group (not considering phylogeny and significant tail-dep.)
# 78 sp. in total (36 LT, 42 UT)
source(here("R/summarize_res.R")) 

#--------- AVONET: traits, bodymass, and phylogeny data for bird species ----------
# file to get traits for all 78 species 
source(here("R/get_birdtraits_from_AVONET.R")) 
source(here("R/get_bodymass_from_AVONET.R")) 
source(here("R/summarize_traits_mass.R"))

source(here("R/get_birdspecies_phylotree.R")) # file to get matched names for all 78 species from BirdTree
#-----------------------------------------------------
# Now consider sp. which shows significant tail-dep spatial synchrony
# based on 75% CI there are 59 out of 78 sp. showing sig results

source(here("R/function_to_testing_abund.R"))
source(here("R/function_to_testsig_tas.R"))

source(here("R/function_to_testsig_pr.R"))
source(here("R/function_to_testsig_pr_avgAprtoAug.R"))

source(here("R/function_to_testsig_tas_avgAprtoAug.R"))
source(here("R/function_to_testsig_tas_avgAprtoJuly.R"))
source(here("R/function_to_testsig_tas_avgMaytoJuly.R"))


source(here("R/distance_sigtaildep_abund_clim.R"))

source(here("R/phylopath_sig75_0-250km.R"))
source(here("R/phylopath_sig75_0-100km.R"))
source(here("R/phylopath_sig75_100-250km.R"))

source(here("R/test_phylogenetic_signal_sig75tree.R")) # maybe put in suppmat

source(here("R/get_migratory_status.R"))
source(here("R/plot_regression.R"))
#=======================





