############################################################
# conceptual figure: copula plot and corresponding timeseries plot
############################################################
#rm(list=ls())
library(VineCopula)
library(copula)
library(here)

source(here("R/ExtremeTailDep.R"))
set.seed(seed=101)

np<-1000000
srho_common<-0.875 #7/8
xtc_r<-retd(n=np,d=2,rl=1,mn=0,sdev=1) # copula with extreme right tail

# make marginals uniform 
xtc_r<-pnorm(xtc_r)
xtc_l<-1-xtc_r

N<-40
time<-1:N

pdf(here("RESULTS/conceptual_fig_cop.pdf"), height=6, width=3)
op<-par(mfrow=c(3,1),mar=c(1,1,1,1),mgp=c(3,1,0),pty="s")

#frank copula
cf<-iRho(frankCopula(),0.875)
f_cop<-BiCopSim(N,family = 5,par=cf)


#=========== pearson correlation ===============
scor_f<-cor.test(f_cop[time,1],f_cop[time,2], method="pearson")
scor_f$estimate

scor_lt<-cor.test(xtc_l[time,1],xtc_l[time,2], method="pearson")
scor_lt$estimate

scor_ut<-cor.test(xtc_r[time,1],xtc_r[time,2], method="pearson")
scor_ut$estimate

#=================================================
# copula plot
# symmetric tail
plot(f_cop[,1],f_cop[,2],col="black",pch=21,
     bg=rgb(0,0,0,0.3),cex=1.5,xlim=c(0,1),ylim=c(0,1),
     xaxt="n",yaxt="n",ann=F,xlab="metapopulation \n abundance, site 1", 
     ylab="metapopulation \n abundance, site 2") # frank

# positive cor with LT
plot(xtc_l[time,1],xtc_l[time,2],col="red",bg=rgb(1,0,0,0.3),pch=21,
     xlim=c(0,1),ylim=c(0,1),
     cex=1.5,
     xaxt="n",yaxt="n",ann=F,xlab="metapopulation \n abundance, site 1", 
     ylab="metapopulation \n abundance, site 2")

# positive cor with UT
plot(xtc_r[time,1],xtc_r[time,2],col="blue",bg=rgb(0,0,1,0.3),pch=21,
     cex=1.5,
     xaxt="n",yaxt="n",ann=F,xlab="metapopulation \n abundance, site 1", 
     ylab="metapopulation \n abundance, site 2")



par(op)
dev.off()


#=========== timeseries plot =================================

pdf(here("RESULTS/conceptual_fig_cop_ts.pdf"), height=6, width=6)
op<-par(mfrow=c(2,2),mar=c(1,1,1,1),mgp=c(1,1,0),pty="s")

# symmetric tail
plot(1:N,f_cop[,1],col="black",pch=1,
     cex=1.5,
     xaxt="n",yaxt="n",ann=F,type="b") # frank
lines(1:N,f_cop[,2],col=rgb(0,0,0,0.5),pch=16,
     cex=1,
     xaxt="n",yaxt="n",ann=F,type="b")

# positive cor with LT
plot(1:N,xtc_l[time,1],col="red",pch=1,
     cex=1.5,
     xaxt="n",yaxt="n",ann=F,type="b") # frank
lines(1:N,xtc_l[time,2],col=rgb(1,0,0,0.3),pch=16,
      cex=1,
      xaxt="n",yaxt="n",ann=F,type="b")

# positive cor with UT
plot(1:N,xtc_r[time,1],col="blue",pch=1,
     cex=1.5,
     xaxt="n",yaxt="n",ann=F,type="b") # frank
lines(1:N,xtc_r[time,2],col=rgb(0,0,1,0.3),pch=16,
      cex=1,
      xaxt="n",yaxt="n",ann=F,type="b")

par(op)
dev.off()





