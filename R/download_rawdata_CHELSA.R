# this script is showing an example for how to download data from CHELSA 
# we only downloaded climate data for 1996-2019 years
library(here)
library(utils)

#=================== download monthly DATA: observed precipitation DATA ======================================
tb<-read.delim(here("DATA/CHELSA_v2/monthly/pr/envidatS3paths.txt"),header = F)
for(i in 1:nrow(tb)){
  readpath<-tb[i,1]
  readpath<-trimws(readpath)
  destfilename<-trimws(basename(readpath))
  year<-as.integer(substr(destfilename,14,17))
#  if(year%in%c(1979:2019)){
    download.file(url=readpath,
                  destfile= paste(here("DATA/CHELSA_v2/monthly/pr"),"/",destfilename,sep=""),mode="wb")
#  }
}

#=================== download monthly DATA: observed mean temperature DATA ======================================
tb<-read.delim(here("DATA/CHELSA_v2/monthly/tas/envidatS3paths.txt"),header = F)
for(i in 1:nrow(tb)){
  readpath<-tb[i,1]
  readpath<-trimws(readpath)
  destfilename<-trimws(basename(readpath))
  year<-as.integer(substr(destfilename,14,17))
#  if(year%in%c(1979:2019)){
    download.file(url=readpath,
                  destfile= paste(here("DATA/CHELSA_v2/monthly/tas"),"/",destfilename,sep=""),mode="wb")
#  }
}

#=================== download monthly DATA: observed max temperature DATA ======================================
tb<-read.delim(here("DATA/CHELSA_v2/monthly/tasmax/envidatS3paths.txt"),header = F)
for(i in 1:nrow(tb)){
  readpath<-tb[i,1]
  readpath<-trimws(readpath)
  destfilename<-trimws(basename(readpath))
  year<-as.integer(substr(destfilename,14,17))
#  if(year%in%c(1979:2019)){
    download.file(url=readpath,
                  destfile= paste(here("DATA/CHELSA_v2/monthly/tasmax"),"/",destfilename,sep=""),mode="wb")
    
#  }
}

#====================================================================================================


