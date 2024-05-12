# run later below code
#==================== plot across species ==============
rm(list=ls())
library(here)
library(tidyverse)
library(ggpubr)
library(gridExtra)
nbin<-4
dfsig75<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail75sig_summary_0-250Km.csv",sep="")))
dfsig95<-read.csv(here(paste("RESULTS/abundance_spatsyn_nbin_",nbin,"_tail95sig_summary_0-250Km.csv",sep="")))


dfsig75.abs<-dfsig75%>%dplyr::select(abs.tot.td.pr5.sig,fpr5.sig,
                                     abs.tot.td.tas5.sig,ftas5.sig,
                                     abs.tot.td.ab.sig, fab.sig)%>%mutate(tailsigCI="75%")

dfsig95.abs<-dfsig95%>%dplyr::select(abs.tot.td.pr5.sig,fpr5.sig,
                                     abs.tot.td.tas5.sig,ftas5.sig,
                                     abs.tot.td.ab.sig, fab.sig)%>%mutate(tailsigCI="95%")

dfsig.abs<-rbind(dfsig75.abs,dfsig95.abs)
dfsig.abs$tailsigCI<-as.factor(dfsig.abs$tailsigCI)
# plot for absolute tail-dep synchrony
gab1<-ggplot(data=dfsig.abs,
           aes(x=abs.tot.td.pr5.sig,y=abs.tot.td.ab.sig,col=tailsigCI), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+ylim(0,5)+
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Precipitation")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=tailsigCI, fill=tailsigCI))+
  theme_bw()+
  stat_cor(aes(col=tailsigCI,
               label = paste(..r.label..,..rr.label.., ..p.label.., 
                             sep = "*`,`~")),
           label.x = c(0.05,0.05), label.y = c(3.85, 4.15), show.legend = F)+
  stat_regline_equation(label.x = c(0.05,0.05), label.y = c(4.45, 4.75),
                        aes(col=tailsigCI),show.legend = F)+
  scale_color_manual(values=c("#E6AB02", "#A6761D"))+
  scale_fill_manual(values=c("#E6AB02", "#A6761D"))+
  theme(legend.position=c(0.8,0.3),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())
gab1

gab2<-ggplot(data=dfsig.abs,
             aes(x=abs.tot.td.tas5.sig,y=abs.tot.td.ab.sig,col=tailsigCI), 
             add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+ylim(0,5)+
  xlab("Total (absolute) tail-dependent spatial synchrony,\n Temperature")+
  ylab("Total (absolute) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=tailsigCI, fill=tailsigCI))+
  theme_bw()+
  stat_cor(aes(col=tailsigCI,
               label = paste(..r.label..,..rr.label.., ..p.label.., 
                             sep = "*`,`~")),
           label.x = c(0.05,0.05), label.y = c(3.85, 4.15), show.legend = F)+
  stat_regline_equation(label.x = c(0.05,0.05), 
                        label.y = c(4.45, 4.75),
                        aes(col=tailsigCI), show.legend = F)+
  scale_color_manual(values=c("#E6AB02", "#A6761D"))+
  scale_fill_manual(values=c("#E6AB02", "#A6761D"))+
  theme(legend.position=c(0.8,0.3),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())
gab2

pdf(here("Results/plot_regression_0-250km_for_abs_taildep.pdf"),width=8, height=4)
grid.arrange(gab1,gab2,nrow=1)
dev.off()

#----------------------
g1<-ggplot(data=dfsig.abs,
           aes(x=fpr5.sig,y=fab.sig,col=tailsigCI), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+ylim(-1,1)+
  xlab("Directional (Net) tail-dependent spatial synchrony,\n Precipitation")+
  ylab("Directional (Net) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=tailsigCI, fill=tailsigCI))+
  theme_bw()+
  stat_cor(aes(col=tailsigCI,
               label = paste(..r.label..,..rr.label.., ..p.label.., 
                             sep = "*`,`~")),
           label.x = c(-1,-1), label.y = c(-0.8, -0.9), show.legend = F)+
  stat_regline_equation(label.x = c(-1,-1), label.y = c(0.85, 0.95),
                        aes(col=tailsigCI),show.legend = F)+
  scale_color_manual(values=c("#E6AB02", "#A6761D"))+
  scale_fill_manual(values=c("#E6AB02", "#A6761D"))+
  theme(legend.position=c(0.8,0.15),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())
g1

g2<-ggplot(data=dfsig.abs,
           aes(x=ftas5.sig,y=fab.sig,col=tailsigCI), 
           add = "reg.line")+
  geom_point(pch=19, alpha=0.5)+
  xlab("Directional (Net) tail-dependent spatial synchrony,\n Temperature")+
  ylab("Directional (Net) tail-dependent spatial synchrony,\n Abundance")+
  geom_smooth(method="lm", se=T,aes(col=tailsigCI, fill=tailsigCI))+
  theme_bw()+
  stat_cor(aes(col=tailsigCI,
               label = paste(..r.label..,..rr.label.., ..p.label.., 
                             sep = "*`,`~")),
           label.x = c(-1,-1), label.y = c(-0.6, -0.75), show.legend = F)+
  stat_regline_equation(label.x = c(-1,-1), label.y = c(-0.25, -0.4),
                        aes(col=tailsigCI),show.legend = F)+
  scale_color_manual(values=c("#E6AB02", "#A6761D"))+
  scale_fill_manual(values=c("#E6AB02", "#A6761D"))+
  theme(legend.position=c(0.8,0.85),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),legend.background=element_blank())
g2

pdf(here("Results/plot_regression_0-250km.pdf"),width=8, height=4)
grid.arrange(g1,g2,nrow=1)
dev.off()

# with phylogeny considered
# pr5_tot abund absolute
fit95<-readRDS(here("RESULTS/model_phylolm_sig95_0-250km/model_tas5_abs_td/model_est_phylolm.RDS"))
fit75<-readRDS(here("RESULTS/model_phylolm_sig75_0-250km/model_tas5_abs_td/model_est_phylolm.RDS"))

plot(fit95, pch=1, col=rgb(red=231/255,green=41/255,blue=138/255,alpha=0.8), 
     ylim=c(-0.1,4), xlim=c(0,4))
par(new=TRUE)
plot(fit75, pch=1, col=rgb(red=117/255,green=112/255,blue=179/255,alpha=0.8), 
     ylim=c(-0.1,4), xlim=c(0,4))


col2rgb("#E7298A")



