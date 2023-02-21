# This plotter function is to visualize a matrix (lower triangular part)
# Input:
#     z : a matrix
#     posnI_ind : a matrix containing the row and col indices of z for which z indicates indep values
#     posnN_ind : a matrix containing the row and col indices of z for which z indicates sig. neg cor values
#     posn_notneeded: a matrix of cell indices to not show in the plot
#     colrange : a vector containing min and max value of the color range
#     tl.cex,cl.cex,line : numeric values for size of text in the labels and colorbar

library(corrplot)

mycorrplot<-function(z,posnI_ind,posnN_ind,colrange,tl.cex,cl.cex){
  
  col1 <- colorRampPalette(c("blue","white","red")) 
  
  z[is.na(z)]<-mean(colrange)
  diag(z)[1]<-colrange[1] # just to ensure that plot always have specific colorbar range even 
  diag(z)[2]<-colrange[2]     # though all entries are either +ve or -ve
  
  colnames(z)[1]<-rownames(z)[1]
  
  corrplot(z,is.corr = F,col=col1(100),method="color",addgrid.col = "black",type="lower",
           diag=F,bg = "white",tl.cex=tl.cex,tl.col = "black",tl.offset = 2.4,tl.pos="ld",
           cl.cex = cl.cex,col.lim = colrange,mar=c(0,0,0,0),cl.length=7,tl.srt=30,
           cl.align.text = "l",cl.ratio = 0.1)
  
  # colorize as black for diagonal indices
  Dg <- matrix(NA,nrow(z),ncol(z))
  diag(Dg)<- 1 
  
  corrplot(Dg, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black",type="lower",
           bg = "transparent", tl.col = "transparent",col="black",method="color")
  
  #colorize as yellow for indep posn indices
  if(dim(posnI_ind)[1]!=0){
    
    I <- matrix(NA,nrow(z),ncol(z))
    I[posnI_ind]<- 1 
    
    corrplot(I, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black",type="lower",
             bg = "transparent", tl.col = "transparent",col="yellow",method="color")
  } 
  
  #colorize as green for -ve correlated (siginificantly) posn indices
  if(dim(posnN_ind)[1]!=0){
    
    N <- matrix(NA,nrow(z),ncol(z))
    N[posnN_ind]<- -1 
    
    corrplot(N, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black",type="lower",
             bg = "transparent", tl.col = "transparent",col="green",pch=16,pch.cex = 2)
  }
  
  # cells with -1 not needed to show
  #if(dim(posn_notneeded)[1]!=0){
  #  znn<-matrix(NA,nrow(z),ncol(z))
  #  znn[posn_notneeded]<- -1
  #  corrplot(znn, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black",type="lower",
  #           bg = "transparent", tl.col = "transparent",col="black",method="color")
  #}
  
}
