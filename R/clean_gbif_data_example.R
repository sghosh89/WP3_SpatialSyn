# This is a script to show an example for a given species how you can clean gbif data 
# you downloaded from gbif

# say, for example, you downloaded gbif occurence data for species "Acanthis flammea" when present
# and with georeferenced info and saved in "DATA/birds_gbif_data/" folder as "raw_Acanthis flammea.csv"

# Below I am showing you how to clean that raw file and save in "DATA/birds_gbif_data/cleaned/" folder

library(here)
library(rgbif)
library(tidyverse)
library(dplyr)
library(countrycode)
library(CoordinateCleaner)
`%notin%` <- Negate(`%in%`)

s<-"Acanthis flammea"
filename<-here(paste("DATA/birds_gbif_data/raw_",s,".csv",sep=""))
x<-read.csv(filename,sep="\t",row.names = NULL)
y <- x%>%mutate(species=s)%>%dplyr::select( c(species,
                                              scientificName,decimalLongitude, decimalLatitude,
                                              basisOfRecord,
                                              issue,
                                              year,month,day,countryCode))%>%
  mutate(occurrenceStatus="PRESENT")%>%rename(issues=issue)

y<-y%>%filter(basisOfRecord%notin%c("FOSSIL_SPECIMEN","PRESERVED_SPECIMEN"))%>%rename(countryCode2=countryCode)
# get 3 letter countrycode
y<-y%>%mutate(countryCode = countrycode(y$countryCode2, origin =  'iso2c', destination = 'iso3c'))


y<-y%>%filter(!is.na(countryCode))
# just to ensure
y$decimalLatitude<-as.numeric(y$decimalLatitude)
y$decimalLongitude<-as.numeric(y$decimalLongitude)
y<-y%>%filter(!is.na(decimalLongitude))%>%filter(!is.na(decimalLatitude))

# automated fixing the issues
flags <- clean_coordinates(x = y, lon = "decimalLongitude", lat = "decimalLatitude", 
                           countries = "countryCode",  species = NULL,
                           tests = c("capitals", "centroids", 
                                     "equal", "gbif", "zeros", 
                                     "countries", "seas"), seas_ref = buffland)

print(summary(flags))
y<-y[flags$.summary,]# cleaned data
write.csv(y,here(paste("DATA/birds_gbif_data/cleaned/",s,".csv",sep="")),row.names = F)

# We followed the above procedure for all bird species occurred in our dataset, but not showing you 
# here as gbif data set for >250 bird species needs a lot of storage space in your machine

