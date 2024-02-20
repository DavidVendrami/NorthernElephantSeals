glc<-"#253741"
fitn<-"#DEC102" #"#FF7F00FF"
rgl<-"#FF7F00FF" #"#E41A1CFF"
dgl<-"#984EA3FF"
igl<-"#4DAF4AFF"
#postc<-"#E41A1CFF"
sims<-"grey60"
sims2<-"grey70"
sims3<-"grey50"
tra<-0.2
tra2<-0.075
rectr<-0.1
mtra<-0.6

highl<-c("#984EA3FF")
#highl<-c("grey")

# The following function was written by someone else. Unfortunately I forgot where i found it
# and I can't attribute it to whom wrote it. If you recognize this as yours, please contact me 
# and I'll be happy to update these lines accordingly. 
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.02,0.9),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, adj=0,...)
}

pdf("Loads_Other_Updated.pdf",width=12, height=8)
par(mfrow=c(3,4))

# All input data here below are output from the various SLiM simulations.

data<-read.table("WEx_Standard/GL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,20),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,20),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,20),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(glc,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="a) total load - exons", location="topleft",cex=1.5)

#############

data<-read.table("WEx_Standard/RGL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(rgl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="b) realized load - exons", location="topleft",cex=1.5)

#############

data<-read.table("WEx_Standard/MGL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(igl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="c) inbreeding load - exons", location="topleft",cex=1.5)

##########

data<-read.table("WEx_Standard/DL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="", axes=F,xaxt='n',xlab="",main="",col=transp(dgl,mtra),lwd=3)

box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="d) drift load - exons", location="topleft",cex=1.5)

#########

data<-read.table("Boot_estimates/GL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,20),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,20),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,20),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(glc,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="e) total load - bootstrap", location="topleft",cex=1.5)

######

data<-read.table("Boot_estimates/RGL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(rgl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="f) realized load - bootstrap", location="topleft",cex=1.5)

#####

data<-read.table("Boot_estimates/MGL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(igl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="g) inbreeding load - bootstrap", location="topleft",cex=1.5)

#####

data<-read.table("Boot_estimates/DL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="", axes=F,xaxt='n',xlab="",main="",col=transp(dgl,mtra),lwd=3)

box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="h) drift load - bootstrap", location="topleft",cex=1.5)


data<-read.table("WGS_Bott10/GL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,20),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,20),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,20),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(glc,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="i) total load - WGS", location="topleft",cex=1.5)

#########

data<-read.table("WGS_Bott10/RGL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(rgl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="j) realized load - WGS", location="topleft",cex=1.5)

#######

data<-read.table("WGS_Bott10/MGL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(igl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="k) inbreeding load - WGS", location="topleft",cex=1.5)

#####

data<-read.table("WGS_Bott10/DL_all.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="", axes=F,xaxt='n',xlab="",main="",col=transp(dgl,mtra),lwd=3)

box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.5, line = 2)
put.fig.letter(label="l) drift load - WGS", location="topleft",cex=1.5)

############################################################

dev.off()
