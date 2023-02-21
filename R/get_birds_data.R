#------------------------------------------------------
#rm(list=ls())
#-----------------------
library(here)
library(dplyr)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
yr_threshold<-20 # we will consider only those routes (i.e., sites) sampled min 20 years
#---------------------
# read the data
xroutes<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/routes.csv")) # route meta data
xweather<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/weather.csv")) # weather data
# BBS provides local weather data, but to keep consistency from different source we are using CHELSA
# it's a problem to read this text file: not in good format
xsplist<-read.delim2(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.txt"),
                     header=F,sep=" ",fileEncoding="latin1")
# read species list files
x<-read.fwf(file = here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.txt"),
            widths = c(6, 
                       7, 
                       50, 51, 50, 
                       19,51,57,65), 
            strip.white = T,fileEncoding="latin1")

colnames(x)<-c("Seq",
               "AOU",
               "English_Common_Name","Freanch_Common_Name","Spanish_Common_Name",
               "ORDER","Family","Genus","Species")
x$ScientificName<-paste(x$Genus,x$Species,sep=" ")

write.csv(x,here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.csv"),row.names = F)

#-----------  I am choosing the species level info for each route (1997-2019) from 50-StopData folder --------------
# for StateNum = 02,Alabama; 03,Alaska; 04,Alberta; 06,Arizona; 07,Arkansas
f1<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty1.csv"))
# for StateNum = 11,BritishColumbia; 14,California; 17,Colorado; 18,Connecticut; 21,Delaware
f2<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty2.csv"))
# for StateNum = 25,Florida;27,Georgia; 33,Idaho; 34,Illinois; 35,Indiana; 36,Iowa
f3<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty3.csv"))
# for StateNum = 38,Kansas; 39,Kentucky;42,Louisiana; 43,Northwest Territories; 44,Maine; 45,Manitoba; 46,Maryland;47,Massachusetts;
f4<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty4.csv"))
# for StateNum = 49,Michigan; 50,Minnesota; 51,Mississippi; 52,Missouri;53,Montana; 54,Nebraska; 55,Nevada; 56,New Brunswick; 57,Newfoundlandand Labrador;
f5<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty5.csv"))
# for StateNum = 58,New Hampshire; 59,New Jersey; 60,New Mexico; 61,NewYork; 63,North Carolina;
f6<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty6.csv"))
# for StateNum = 62,Nunavut; 64,North Dakota; 65,Nova Scotia;66,Ohio; 67,Oklahoma; 68,Ontario;
f7<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty7.csv"))
# for StateNum = 69,Oregon; 72,Pennsylvania; 75,PrinceEdward Island; 76,Quebec; 77,Rhode Island; 79,Saskatchewan; 80,SouthCarolina;
f8<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty8.csv"))
# for StateNum = 81,South Dakota; 82,Tennessee; 83,Texas; 85,Utah; 87,Vermont
f9<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty9.csv"))
# for StateNum = 88,Virginia; 89,Washington; 90,West Virginia; 91,Wisconsin; 92,Wyoming;93,Yukon
f10<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty10.csv"))

ffull<-rbind(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)
saveRDS(ffull,here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty1to10.RDS"))
#ffull<-readRDS(here(DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/50-StopData/1997ToPresent_SurveyWide/fifty1to10.RDS"))

#------------------------ now data cleaning --------------------------------
#Note that the dataset is complete only for the period from 1997 to present; most earlier data
#(1966-1996) are only available in the 10-stop summary files (in the "States" folder, online data).
fshort<-ffull%>%filter(Year%in%c(1997:2019)) # filter for 1997-2019
colnames(fshort)
sum1to50<-apply(fshort[,8:57],MARGIN=1,FUN=sum)
fshort<-fshort%>%
  mutate(Stop1to50=sum1to50)%>% # total count for all 50 stop on each route
  dplyr::select(RouteDataID,CountryNum,StateNum,Route,RPID,Year,AOU,Stop1to50) # selecting some columns

# now, consider the RouteDataID which were consistently maintained BBS protocol
xweather_good<-xweather%>%filter(Year%in%c(1997:2019))%>%
  filter(RunType==1)
fshort<-fshort%>%filter(RouteDataID%in%xweather_good$RouteDataID)

saveRDS(fshort,here("DATA/for_BBS/wrangled_data/data1997to2019_consistentprotocol.RDS"))
zz<-head(fshort,10)

fshort$uRID<-paste(fshort$CountryNum,"_",fshort$StateNum,"_",fshort$Route,sep="")

#How many routes are there?
length(unique(fshort$uRID))#4687

# one given route can be among many states
# 1227 routes are sampled atleast for 20 years, we will consider only those
c2<-fshort%>%group_by(uRID)%>%summarise(nyr=n_distinct(Year))%>%ungroup()%>%filter(nyr>=yr_threshold)

xroutes<-xroutes%>%dplyr::select(CountryNum, StateNum, Route, RouteName, Latitude, Longitude, Stratum)
xroutes<-xroutes%>%mutate(uRID=paste(CountryNum, StateNum, Route,sep="_"))

route_20yrs_metadata<-left_join(c2,xroutes,by="uRID")
saveRDS(route_20yrs_metadata,here("DATA/for_BBS/wrangled_data/route_20yrs_metadata.RDS"))

# Now we need to compute the pairwise geographic distance of this 1227 sites
library(geosphere)
dist<-distm(x=route_20yrs_metadata[,c("Longitude","Latitude")], fun=distGeo)
dist<-dist/1000 # distance in Km
dim(dist)
rownames(dist)<-colnames(dist)<-route_20yrs_metadata$uRID
# dist is a symmetric matrix, with pairwise distance between sites, in Km
range(dist)
dist[upper.tri(dist, diag=T)]<-NA
saveRDS(dist,here("DATA/for_BBS/wrangled_data/pairwise_distance_km.RDS"))
# some visualization
hist(dist,100)
range(dist,na.rm=T)
# we will choose some distance category later: <250km, 250-500km,
# 500-1000 km, 1000-2000 km, 2000-3000km, 3000-7000km

# Now we also need to get a species presence-absence matrix across sites
# i.e., for any given species, get a matrix where sites along columns,
# years across rows

fshort<-readRDS(here("DATA/for_BBS/wrangled_data/data1997to2019_consistentprotocol.RDS"))
fshort$uRID<-paste(fshort$CountryNum,"_",fshort$StateNum,"_",fshort$Route,sep="")
fshort_20yr<-fshort%>%filter(uRID%in%route_20yrs_metadata$uRID)
saveRDS(fshort_20yr,here("DATA/for_BBS/wrangled_data/data1997to2019_consistentprotocol_min20yr.RDS"))

nsp<-length(unique(fshort_20yr$AOU))# 652 species recorded
nyr<-length(unique(fshort_20yr$Year)) # 23 years sampled
range(fshort_20yr$Year)# 1997-2019
length(unique(fshort$uRID))
nsite<-length(unique(fshort_20yr$uRID)) #1227 sites
usp<-sort(unique(fshort_20yr$AOU)) # 652 unique species
  
abund_array<-array(data=0, dim=c(nyr,nsite,nsp))
dimnames(abund_array)[[1]]<- as.character(sort(unique(fshort_20yr$Year)))
dimnames(abund_array)[[2]]<- route_20yrs_metadata$uRID
dimnames(abund_array)[[3]]<- as.character(usp)

species_absentinfo<-data.frame(AOU=NA*numeric(nsp),nAbsentSitesAllyr=NA*numeric(nsp))

# now fill in the array

for(i in 1:nsp){
  
  sp<-usp[i]
  zz<-fshort_20yr%>%filter(AOU==sp)%>%dplyr::select(uRID,Year,Stop1to50)
  notobsyr<-setdiff(c(1997:2019),zz$Year)
  zz_ref<-data.frame(Year=notobsyr)
  newtab<-full_join(zz,zz_ref,by="Year")
  zz0<-newtab %>% complete(uRID, Year, fill= list(Stop1to50=0))
  zz0<-na.omit(zz0)
  z0<-zz0%>%group_split(uRID)
  sitename<-unlist(sapply(z0, function(x) x[1, 1]))
  matchsiteid<-match(sitename,dimnames(abund_array)[[2]])
  myabund<-sapply(z0, function(x) x[, 3])
  
  for(j in 1:length(myabund)){
    siteid<-matchsiteid[j]
    abund_array[,siteid,i]<-myabund[[j]]
    cat(paste("i,j=",i,j,"\n"))
  }
  
  mat<-abund_array[,,i]
  all0<-apply(mat,MARGIN=2,FUN=sum)
  id_all0<-which(all0==0)
  species_absentinfo$AOU[i]<-sp
  species_absentinfo$nAbsentSitesAllyr[i]<-length(id_all0)# number of sites where given species never observed throughout the years
}
dim(abund_array)# year by site by species

saveRDS(abund_array,here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_array.RDS"))


write.csv(species_absentinfo,here("DATA/for_BBS/wrangled_data/data1997to2019_abundance_species_absentinfo.csv"), row.names = F)

xroutes<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/routes.csv")) # route meta data
xroutes$uRID<-paste(xroutes$CountryNum,"_",xroutes$StateNum,"_",xroutes$Route,sep="")
xroutes<-xroutes%>%dplyr::select(uRID,Longitude,Latitude,Stratum)
write.csv(xroutes,here("DATA/for_BBS/wrangled_data/uRID_lonlat_stratum.csv"), row.names = F)

# so here 652 species are filtered because they were at least sampled once in the study period





