# run later below code
#==================== plot across species ==============
rm(list=ls())
library(here)
library(tidyverse)
library(ggpubr)
library(gridExtra)
nbin<-4
dfsig<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail75sig_summary_0-250Km.csv",sep="")))

# plot for absolute tail-dep synchrony
gab1<-ggplot(data=dfsig,
           aes(x=abs.tot.td.pr5.sig,y=abs.tot.td.ab.sig), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.3)+#ylim(c(0,4))+xlim(c(0,4))+
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Precipitation")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", col="black",se=T)+
  theme_bw()+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = 0.85, label.y = 0.75,col="black")+
  stat_regline_equation(label.x = 0.85, label.y = 0.3)+
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
gab1

gab2<-ggplot(data=dfsig,
             aes(x=abs.tot.td.tas5.sig,y=abs.tot.td.ab.sig), 
             add = "reg.line")+
  geom_point(pch=19, alpha=0.3)+#ylim(c(0,4))+xlim(c(0,4))+
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Temperature")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", col="black",se=T)+
  theme_bw()+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = 0.7, label.y = 0.55,col="black")+
  stat_regline_equation(label.x = 0.7, label.y = 0.1)+
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
gab2

pdf(here("Results/plot_regression_0-250km_for_abs_taildep.pdf"),width=8, height=4)
grid.arrange(gab1,gab2,nrow=1)
dev.off()

#----------------------
g1<-ggplot(data=dfsig,
           aes(x=fpr5.sig,y=fab.sig), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.3)+ylim(-1,1.6)+
  xlab("Tail-dependent spatial synchrony,\n Precipitation")+
  ylab("Tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", col="black",se=T)+
  theme_bw()+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = -0.95, label.y = -0.75,col="black")+
  stat_regline_equation(label.x = -0.95, label.y = -0.95)+
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g1

g2<-ggplot(data=dfsig,aes(x=ftas5.sig,y=fab.sig), add = "reg.line")+
  geom_point(pch=19, alpha=0.3)+ylim(-1,1.6)+
  xlab("Tail-dependent spatial synchrony,\n Temperature")+
  ylab("Tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", col="black",se=T)+
  theme_bw()+
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x = -0.95, label.y = -0.55,col="black")+
  stat_regline_equation(label.x = -0.95, label.y = -0.8)+
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g2

pdf(here("Results/plot_regression_0-250km.pdf"),width=7, height=3)
grid.arrange(g1,g2,nrow=1)
dev.off()


