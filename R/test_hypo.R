# We will test some hypotheses related to trait-difference 
# among two groups showing LT and UT spatial synchrony

library(tidyverse)
library(dplyr)
library(here)
library(ggpubr)
library(gridExtra)
library(PupillometryR)

chosen_rad<-c(0,250)

df<-read.csv(here(paste("RESULTS/df_abund_climate_spatsyn_",chosen_rad[1],"_",
chosen_rad[2],"km.csv", sep="")))
df_spmeta<-read.csv(here("RESULTS/species_dietcat_edited.csv"))
df_spmeta<-df_spmeta%>%dplyr::select(AOU, English_Common_Name, ScientificName,PassNonPass)

df<-left_join(df,df_spmeta, by="AOU") # don't rename this dataframe


#=========================  BODY SIZE RELATED TRAITS  ==========================

# first with Bird AvoNet database
dftrait<-read.csv(here("DATA/AVONET/bird_traits_from_AVONET.csv"))
dftrait<-dftrait%>%dplyr::select(spname, possible_sp, 
                                 Beak.Length_Culmen, Beak.Width, Beak.Depth,
                                 Tarsus.Length, 
                                 Wing.Length, Kipps.Distance, Hand.wing.Index,
                                 Tail.Length)
dftrait<-dftrait%>%group_by(spname,possible_sp)%>%
  summarize(meanBeak.LengthCulmen=mean(Beak.Length_Culmen,na.rm=T),
            meanBeak.Width=mean(Beak.Width,na.rm=T),
            meanBeak.Depth=mean(Beak.Depth,na.rm=T),
            meanTarsus.Length=mean(Tarsus.Length,na.rm=T),
            meanWing.Length=mean(Wing.Length,na.rm=T),
            meanKipps.Distance=mean(Kipps.Distance,na.rm=T),
            meanHWI=mean(Hand.wing.Index,na.rm=T),
            meanTail.Length=mean(Tail.Length,na.rm=T),)%>%ungroup()

dft1<-dftrait%>%filter(spname==possible_sp)%>%dplyr::select(-spname)
dft2<-dftrait%>%filter(spname!=possible_sp)%>%dplyr::select(-possible_sp)

df2<-left_join(df,dft1,by=c("ScientificName"="possible_sp"))
df2<-left_join(df2,dft2,by=c("ScientificName"="spname"))

df2$meanBeak.LengthCulmen<-coalesce(df2$meanBeak.LengthCulmen.x,df2$meanBeak.LengthCulmen.y)
df2$meanBeak.Width<-coalesce(df2$meanBeak.Width.x,df2$meanBeak.Width.y)
df2$meanBeak.Depth<-coalesce(df2$meanBeak.Depth.x,df2$meanBeak.Depth.y)
df2$meanTarsus.Length<-coalesce(df2$meanTarsus.Length.x,df2$meanTarsus.Length.y)
df2$meanWing.Length<-coalesce(df2$meanWing.Length.x,df2$meanWing.Length.y)
df2$meanKipps.Distance<-coalesce(df2$meanKipps.Distance.x,df2$meanKipps.Distance.y)
df2$meanHWI<-coalesce(df2$meanHWI.x,df2$meanHWI.y)
df2$meanTail.Length<-coalesce(df2$meanTail.Length.x,df2$meanTail.Length.y)

# data with 8 traits from AvoNet
mydftrait<-df2%>%dplyr::select(AOU, English_Common_Name, ScientificName, 
                               fLU_ab,fLU_pr, fLU_tas, fLU_tasmax, fLU_tasmin,Diet.5Cat, 
                               IUCN_status, tail, PassNonPass,
                               meanBeak.LengthCulmen,meanBeak.Width,meanBeak.Depth,
                               meanTarsus.Length,meanWing.Length,
                               meanKipps.Distance,meanHWI,
                               meanTail.Length)

plot_type_bodytrait<-function(mydftrait,traits){
  
  gplist<-vector(mode = "list", length = length(traits))
  
  for(i in 1:length(traits)){
    
      id<-which(colnames(mydftrait)%in%c("tail","Diet.5Cat",
                                         "fLU_ab","fLU_pr",
                                         "fLU_tas","fLU_tasmax","fLU_tasmin",
                                         traits[i]))
      tempo<-mydftrait[id]
      tempo<-rename(tempo,DT=Diet.5Cat,y=traits[i],x=fLU_ab)
      gp<-ggplot(tempo, aes(x=x, y=y, col=tail, shape=tail)) + 
        geom_point(alpha=0.5)+geom_smooth(se=F,method="lm")+
        labs(title="",x="spat syn (abund)",y=traits[i])+
        theme_bw()+theme(legend.position = "none", legend.box = "horizontal", 
                         axis.text.x = element_text(size = 15),
                         axis.text.y = element_text(size = 15),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank())+
        stat_cor(method = "spearman") 
        #stat_cor(aes(label = paste(after_stat(p.label), sep = "~")))
   
    
    gplist[[i]]<-gp
    
  }
  
  return(gplist)
}

traits<-c("meanBeak.LengthCulmen","meanBeak.Width","meanBeak.Depth",
          "meanTarsus.Length","meanWing.Length","meanKipps.Distance",
          "meanHWI","meanTail.Length")

# but here we can see in extreme conditions because of some environmental filtering species with better traits exist
res<-plot_type_bodytrait(mydftrait=mydftrait,traits=traits)
grid.arrange(res[[1]], res[[2]], res[[3]], res[[4]], res[[5]], res[[6]], res[[7]], res[[8]],  
             layout_matrix = rbind(c(1, 2, 3, 4),
                                   c(5, 6, 7, 8)), nrow=2)

pdf(here(paste("RESULTS/species_bodytraits_vs_spatsynabundvalue_",
               chosen_rad[1],"_",chosen_rad[2],"km.pdf",sep="")), width = 12, height = 6)
grid.arrange(res[[1]], res[[2]], res[[3]], res[[4]], res[[5]], res[[6]], res[[7]], res[[8]],  
             layout_matrix = rbind(c(1, 2, 3, 4),
                                   c(5, 6, 7, 8)), nrow=2)
dev.off()
#=========================  BODY MASS   ==========================

dfmass<-read.csv(here("DATA/AVONET/bird_bodymass_from_AVONET.csv"))# 13614 species
dfm<-left_join(df,dfmass,by=c("ScientificName"="Species"))
dfm1<-dfm%>%filter(!is.na(mass_gm)) # 251 sp.
dfm2<-dfm%>%filter(is.na(mass_gm)) # 11 sp. - need to know synonym

dftrait<-read.csv(here("DATA/AVONET/bird_traits_from_AVONET.csv"))
dftrait<-dftrait%>%filter(spname!=possible_sp)%>%distinct(spname,.keep_all = T)
dftrait<-dftrait%>%dplyr::select(spname,possible_sp)
dfm2<-left_join(dfm2,dftrait,by=c("ScientificName"="spname"))
id<-which(dfm2$ScientificName%in%c("Setophaga coronata coronata","Setophaga coronata audoboni"))
dfm2$possible_sp[id]<-"Setophaga coronata" # for these subspecies mass info not available, so putting species level info
dfm2<-left_join(dfm2,dfmass,by=c("possible_sp"="Species"))
dfm2$mass_gm<-coalesce(dfm2$mass_gm.x,dfm2$mass_gm.y)
dfm2$source<-coalesce(dfm2$source.x,dfm2$source.y)
dfm2<-dfm2%>%dplyr::select(-source.x,-source.y,-mass_gm.x,-mass_gm.y)

dfm1$possible_sp<-dfm1$ScientificName
dfm2<-dfm2%>%dplyr::select(colnames(dfm1))

# combine
dfm_birds<-rbind(dfm1,dfm2)%>%arrange(AOU)
dfm_birds2<-dfm_birds%>%dplyr::select(ScientificName, possible_sp, mass_gm, source)
mydftrait_w_mass<-left_join(mydftrait,dfm_birds2,by=c("ScientificName"="ScientificName"))

# save file
write.csv(mydftrait_w_mass,here("RESULTS/df_abund_climate_spatsyn_0_250km_with_speciestraits_mass.csv"),row.names=F)

res<-plot_type_bodytrait(mydftrait=dfm_birds,traits="mass_gm")
print(res[[1]])

pdf(here(paste("RESULTS/species_bodymass_vs_spatsynabundvalue_",
               chosen_rad[1],"_",chosen_rad[2],"km.pdf",sep="")), width = 3, height = 3)
print(res[[1]])
dev.off()












