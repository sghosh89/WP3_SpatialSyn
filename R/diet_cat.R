library(here)
library(tidyverse)
`%notin%` <- Negate(`%in%`)

x<-read.csv(here("DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/SpeciesList_edited.csv"))                                       
y<-read.delim(here("DATA/for_BBS/raw_data/EltonTraits/BirdFuncDat.txt"))

getcat<-function(x){
  res<-which(x>50)
  if(length(res)==0){
    resy<-"mixed"
  }else{
    resy<-names(x)[res]
    resy<-str_split(resy,pattern="\\.",n=2)[[1]][2]
  }
  return(resy)
}
z1<-y%>%select(ForStrat.watbelowsurf,
               ForStrat.wataroundsurf,
               ForStrat.ground,ForStrat.understory,
               ForStrat.midhigh,ForStrat.canopy,
               ForStrat.aerial)
z11<-apply(z1,MARGIN=1,FUN=sum)
table(z11)# not all 100 as sometimes rounded number would give 99 

s<-apply(z1,MARGIN = 1, FUN=getcat)
y$Strat.7Cat<-s
y<-y%>%select(SpecID,PassNonPass,IOCOrder,BLFamilyLatin,BLFamilyEnglish,BLFamSequID,Taxo,
              Scientific,English,
              Diet.5Cat,Diet.Certainty,
              Strat.7Cat,ForStrat.SpecLevel,
              Nocturnal,PelagicSpecialist)

# matching by scientific name: 534 out of 756 matched
z<-inner_join(x,y, by=c("ScientificName"="Scientific"))
table(z$Diet.5Cat) # 5 diet categories
# PlantSeed: Plant and Seeds, 
# FruiNect: Fruits and Nectar, 
# Invertebrate: Invertebrates, 
# VertFishScav: Vertebrates and Fish and Carrion, 
# Omnivore: score of <= 50 in all four categories
table(z$Strat.7Cat)
# 7 categories based on prevalence of Foraging
# ForStrat-watbelowsurf: Prevalence of Foraging below the water surfaces
# ForStrat-wataroundsurf: Prevalence of Foraging on or just (<5 inches) below water surface
# ForStrat-ground: Prevalence of Foraging on ground
# ForStrat-understory: Prevalence of Foraging below 2m in understory in forest, forest edges, bushes or shrubs
# ForStrat-midhigh: Prevalence of Foraging in mid to high levels in trees or high bushes (2m upward), but below canopy
# ForStrat-canopy: Prevalence of Foraging in or just above (from) tree canopy
# ForStrat-aerial: Prevalence of Foraging well above vegetation or any structures

z0<-anti_join(x,y, by=c("ScientificName"="Scientific"))# 222 out of 756 not matched

# 71% matched after scientific name matching 

##########################################################
# Now, we will check if some unmatched (no scientific name matching) entries could be matched based on english name?
ynoScientificmatch<-anti_join(y,z, by=c("Scientific"="ScientificName"))
zEngName<-inner_join(z0,ynoScientificmatch, by=c("English_Common_Name"="English")) # 117 matched


# species found till now
z<-rename(z, mixed_name=English)
zEngName<-rename(zEngName, mixed_name= Scientific)

dftemp<-rbind(z,zEngName)

# check if the 373 species AOU you choose here for spat syn calculation are within the z dataframe or not?
df<-read.csv(here("DATA/for_BBS/wrangled_data/data1979to2019_abundance_species_w_minimum2sites.csv"))

dfall<-left_join(df,dftemp, by="AOU") # 9 species data not found, we will fill in manually
df0<-dfall[which(is.na(dfall$Seq)),]
write.csv(dfall, here("RESULTS/species_dietcat.csv"),row.names = F)

#df0<-dfall[which(is.na(dfall$Seq)),]

# I filled in the edited csv mainly for diet category and for those 9 species 
# from df0 dataframe, x

#x<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
#old<-read.csv(here("RESULTS/old.csv"))
#old<-old%>%dplyr::select(AOU, IUCN_status)
#y<-left_join(x,old,by="AOU")
#y<-y%>%dplyr::select(-IUCN.status)
#write.csv(y, here("RESULTS/species_dietcat_edited.csv"),row.names = F)

