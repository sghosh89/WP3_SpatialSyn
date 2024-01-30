# Here we will extract data for birds species count before 1997
# this data are avialable in states folder as 10 stop summary
xweather<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/weather.csv")) # weather data

resloc<-here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/States")

filelist <- list.files(path = resloc, pattern='.csv$', all.files=TRUE, full.names=FALSE)

readfile<-paste(resloc,filelist,sep="/")
fl1<-read.csv(readfile[1]);fl2<-read.csv(readfile[21])
all(colnames(fl1)==colnames(fl2))==T # that means all csv files in same format

files<-c()
for(i in 1:62){
  tempo<-read.csv(readfile[i])
  tempo<-tempo%>%filter(Year%in%c(1979:1996))
  files<-rbind(files,tempo)
}
files<-files%>%dplyr::select(RouteDataID,CountryNum,StateNum,Route,RPID,Year,AOU,
                             Stop1to50=SpeciesTotal)

xweather_good<-xweather%>%filter(Year%in%c(1979:1996))%>%
  filter(RunType==1)
files<-files%>%filter(RouteDataID%in%xweather_good$RouteDataID)
saveRDS(files,here("DATA/for_BBS/wrangled_data/data1979to1996_consistentprotocol.RDS"))


