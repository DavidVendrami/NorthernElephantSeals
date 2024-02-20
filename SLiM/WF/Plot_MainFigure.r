# Plot

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
siz<-"#1f78b4"

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

nes<-read.table("Ne_RAD_table.txt", h=T)

pdf("WFsimulations_Ne_V3.pdf",width=7, height=7.5)
par(mfrow=c(2,2))

######

par(mar = c(5.1, 4.1, 4.1, 3.1))
data<-read.table("GL_all.txt",h=F) # Simulation output

plot(0:24,nes$N,type='l',col="white",yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",lwd=2)
polygon(c(0:24,24:0),c(rep(0,25),nes$N[25:1]),col=transp(siz,.3),border=transp(siz,.3))
#plot(0:24,nes$N,type='l',col=transp(siz,.3),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",lwd=2)
axis(4,c(-10000,0,5000,10000,15000),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
#axis(4,c(-10000,0,5000,10000,15000),labels=c(as.character(c(-10000,0,5,10,15))),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
par(xpd=T)
text(29,6968.085,expression(italic('N')[e]),srt=270,cex=1.2)
#text(29,6968.085,expression(italic('N')[e]*" (x 10"^3*")"),srt=90,cex=1.2)

par(xpd=F)
par(new=T)

ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,20),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.2)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,20),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,20),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(glc,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.2, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.2, line = 2)
put.fig.letter(label="a) total load", location="topleft",cex=1.2)

#######

par(mar = c(5.1, 4.1, 4.1, 3.1))
data<-read.table("RGL_all.txt",h=F) # Simulation output

plot(0:24,nes$N,type='l',col="white",yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",lwd=2)
polygon(c(0:24,24:0),c(rep(0,25),nes$N[25:1]),col=transp(siz,.3),border=transp(siz,.3))
#plot(0:24,nes$N,type='l',col=transp(siz,.3),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",lwd=2)
axis(4,c(-10000,0,5000,10000,15000),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
#axis(4,c(-10000,0,5000,10000,15000),labels=c(as.character(c(-10000,0,5,10,15))),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
par(xpd=T)
text(29,6968.085,expression(italic('N')[e]),srt=270,cex=1.2)
#text(29,6968.085,expression(italic('N')[e]*" (x 10"^3*")"),srt=90,cex=1.2)
par(xpd=F)
par(new=T)

ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.2)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(rgl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.2, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.2, line = 2)
put.fig.letter(label="b) realized load", location="topleft",cex=1.2)

#######

par(mar = c(5.1, 4.1, 4.1, 3.1))
data<-read.table("MGL_all.txt",h=F) # Simulation output

plot(0:24,nes$N,type='l',col="white",yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",lwd=2)
polygon(c(0:24,24:0),c(rep(0,25),nes$N[25:1]),col=transp(siz,.3),border=transp(siz,.3))
#plot(0:24,nes$N,type='l',col=transp(siz,.3),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",lwd=2)
axis(4,c(-10000,0,5000,10000,15000),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
#axis(4,c(-10000,0,5000,10000,15000),labels=c(as.character(c(-10000,0,5,10,15))),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
par(xpd=T)
text(29,6968.085,expression(italic('N')[e]),srt=270,cex=1.2)
#text(29,6968.085,expression(italic('N')[e]*" (x 10"^3*")"),srt=90,cex=1.2)
par(xpd=F)
par(new=T)

ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.2)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(igl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.2, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.2, line = 2)
put.fig.letter(label="c) inbreeding load", location="topleft",cex=1.2)

########

par(mar = c(5.1, 4.1, 4.1, 3.1))
data<-read.table("DL_all.txt",h=F) # Simulation output

plot(0:24,nes$N,type='l',col="white",yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",lwd=2)
polygon(c(0:24,24:0),c(rep(0,25),nes$N[25:1]),col=transp(siz,.3),border=transp(siz,.3))
#plot(0:24,nes$N,type='l',col=transp(siz,.3),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",lwd=2)
axis(4,c(-10000,0,5000,10000,15000),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
#axis(4,c(-10000,0,5000,10000,15000),labels=c(as.character(c(-10000,0,5,10,15))),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
par(xpd=T)
text(29,6968.085,expression(italic('N')[e]),srt=270,cex=1.2)
#text(29,6968.085,expression(italic('N')[e]*" (x 10"^3*")"),srt=90,cex=1.2)
par(xpd=F)
par(new=T)

ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,-(data$V1),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.2)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,-(data[,i]),type='l',ylim=c(0,15),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot(-(ms),type='l',ylim=c(0,15),yaxt='n',ylab="", axes=F,xaxt='n',xlab="",main="",col=transp(dgl,mtra),lwd=3)

box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,5,10,15,20),cex.axis=1.2,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.2, line = 2)
title(ylab = "Lethal equivalents", cex.lab = 1.2, line = 2)
put.fig.letter(label="d) drift load", location="topleft",cex=1.2)

#######

dev.off()

