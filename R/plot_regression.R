# run later below code
#==================== plot across species ==============
library(ggpubr)
nbin<-4
dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail75sig_summary_0-250Km.csv",sep="")))

df$tail75<-factor(df$tail75)
ftmethod<-"lm"

g1<-ggplot(data=df,aes(x=fpr5.sig,y=fab.sig), add = "reg.line")+
  geom_point(pch=21, col="white")+
  #geom_smooth(method=ftmethod,linetype="dashed")+
  xlab("Tail-dep. spatial synchrony, Precipitation")+
  ylab("Tail-dep. spatial synchrony, Abundance")+
  #facet_wrap(~Diet.5Cat)+
  geom_smooth(method=ftmethod, col="black",se=T)+
  theme_bw()+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = -0.7, label.y = 0.8,col="black")+
  stat_regline_equation(label.x = -0.7, label.y = 0.9)+
  geom_point(data=df, aes(x=fpr5.sig,y=fab.sig, col=tail75), alpha=0.3)+
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g1

g2<-ggplot(data=df,aes(x=ftas5.sig,y=fab.sig), add = "reg.line")+
  geom_point(pch=21, col="white")+
  #geom_smooth(method=ftmethod,linetype="dashed")+
  xlab("Tail-dep. spatial synchrony, Temperature")+
  ylab("Tail-dep. spatial synchrony, Abundance")+
  #facet_wrap(~Diet.5Cat)+
  geom_smooth(method=ftmethod, col="black",se=T)+
  theme_bw()+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = -0.7, label.y = 0.8,col="black")+
  stat_regline_equation(label.x = -0.7, label.y = 0.9)+
  geom_point(data=df, aes(x=ftas5.sig,y=fab.sig, col=tail75), alpha=0.3)+
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g2

#=============== plot groupwise: rare on top row, common on bottom row ===========
ftmethod<-"lm"

# ======== NOW WITH SPECIES WHICH ARE SYNCHRONOUSLY RARE (n=138) ===================
# LT synchrony for abundance caused by high pr (+ve cor) and high temp across sites (-ve cor)

df1<-df%>%filter(tail75=="LT")
gLT<-ggplot(data=df1,aes(x=ftas5.sig,y=fab.sig), add = "reg.line")+geom_point(alpha=0.3, col = "red")+
  geom_smooth(method=ftmethod,aes(col="red"))+
  xlab("Tail-dep. spatial synchrony, T")+
  ylab("Tail-dep. spatial synchrony, Abundance")+
  #facet_wrap(~Diet.5Cat)+
  theme_bw()+ 
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = -0.7, label.y = 0.8,col="black")+
  stat_regline_equation(label.x = -0.7, label.y = 0.9)
gLT

# === NOW WITH SPECIES WHICH ARE SYNCHRONOUSLY COMMON ===================
# UT synchrony for abundance caused by high temp across sites (+ve cor)
df2<-df%>%filter(tail75=="UT")
gUT<-ggplot(data=df2,aes(x=ftas5.sig,y=fab.sig))+geom_point(alpha=0.3, col = "darkturquoise")+
  geom_smooth(method=ftmethod, col = "darkturquoise")+
  xlab("Tail-dep. spatial synchrony, T")+
  ylab("Tail-dep. spatial synchrony, Abundance")+
  #facet_wrap(~Diet.5Cat)+
  theme_bw()+ 
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = 0.5, label.y = -0.9,col="black")+
  stat_regline_equation(label.x = 0.5, label.y = -0.8)
gUT

grid.arrange(gLT, gUT, ncol=2)

