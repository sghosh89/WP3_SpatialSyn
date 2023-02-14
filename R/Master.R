#-------------------
path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
#--------------------
source("./get_birds_data.R") # wrangle raw data
source("./prepare_abund_data.R") # for a given species, get abundance data in required format