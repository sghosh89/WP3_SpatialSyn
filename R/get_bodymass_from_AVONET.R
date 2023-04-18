#====================================================================
# This file will get the body mass data for the bird species
# considered for spatial synchrony (0-250Km) results, i.e., 262 species#
#====================================================================

#rm(list=ls())
`%notin%` <- Negate(`%in%`)
library(tidyverse)
library(dplyr)
library(readxl)
library(here)
# first for all bird species

avodat_meta<-read_excel(here("DATA/AVONET/AVONET Supplementary dataset 1.xlsx"),sheet=1)
avomass_BirdLife<-read_excel(here("DATA/AVONET/AVONET Supplementary dataset 1.xlsx"),sheet=2)
avomass_BirdLife<-avomass_BirdLife%>%dplyr::select(Species=Species1,Mass,Mass.Source)
avomass_BirdLife$source<-"BirdLife"

avomass_eBird<-read_excel(here("DATA/AVONET/AVONET Supplementary dataset 1.xlsx"),sheet=3)
avomass_eBird<-avomass_eBird%>%dplyr::select(Species=Species2,Mass,Mass.Source)
avomass_eBird$source<-"eBird"

avomass_BirdTree<-read_excel(here("DATA/AVONET/AVONET Supplementary dataset 1.xlsx"),sheet=4)
avomass_BirdTree<-avomass_BirdTree%>%dplyr::select(Species=Species3,Mass,Mass.Source)
avomass_BirdTree$source<-"BirdTree"

# Combine all data
avomass<-rbind(avomass_BirdLife,avomass_eBird,avomass_BirdTree)
avomass<-avomass%>%arrange(Species)

avomass<-avomass%>%group_by(Species)%>%summarise(mass_gm=mean(Mass,na.rm=T))%>%ungroup()
avomass$source<-"combined: BirdLife, eBird, BirdTree"
write.csv(avomass,here("DATA/AVONET/bird_bodymass_from_AVONET.csv"),row.names=F)
#========================















