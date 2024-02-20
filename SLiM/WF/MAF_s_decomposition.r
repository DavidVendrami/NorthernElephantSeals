glc<-"#253741"
fitn<-"#DEC102" #"#FF7F00FF"
rgl<-"#FF7F00FF" #"#E41A1CFF"
dgl<-"#984EA3FF"
igl<-"#4DAF4AFF"
postc<-"#7198AD"
sims<-"grey60"
sims2<-"grey70"
sims3<-"grey50"
tra<-0.05
tra2<-0.075
rectr<-0.1
mtra<-0.4

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
                     topleft = c(0.02,0.95),
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

pdf("Mutations_trajectories.pdf",width=8.83, height=6.5)
par(mfrow=c(2,3))

################################

# All input data here below are output from the various SLiM simulations.

data<-read.table("Pre_maf05.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,200000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,200000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,200000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(glc,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,100000,200000),labels=c("0","100000","200000"),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Counts", cex.lab = 1.5, line = 2)
put.fig.letter(label="a) rare mutations", location="topleft",cex=1.5)
legend(10,200000,legend=c(expression("AF "<=" 0.05")),
						col = transp(c(glc),.8),lwd=3,lty=1,box.lty=0,cex=1.0,pt.cex=2)

#################################

data<-read.table("Pre_maf10.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(fitn,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(fitn,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,5000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(fitn,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_maf25.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(rgl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(rgl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,5000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(rgl,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_maf50.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(igl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(igl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,5000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(igl,mtra),lwd=3)


par(new=T)

data<-read.table("Pre_maf1.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(postc,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(postc,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,5000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(postc,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_mafFixed.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(dgl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,5000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(dgl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,5000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(dgl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,2500,5000),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Counts", cex.lab = 1.5, line = 2)
put.fig.letter(label="b) more common mutations", location="topleft",cex=1.5)
legend(10,5000,legend=c(expression("0.05 < AF "<=" 0.1"),expression("0.1 < AF "<=" 0.25"),
						expression("0.25 < AF "<=" 0.5"),expression("0.5 < AF "<=" 1"),"AF = 1"),
						col = transp(c(fitn,rgl,igl,postc,dgl),.8),lwd=3,lty=1,box.lty=0,cex=1.0,pt.cex=2)

#############################################################

tots<-read.table("Total_prebott.txt")
data<-read.table("Pre_maf05.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(glc,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_maf10.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(fitn,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(fitn,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(fitn,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_maf25.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(rgl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(rgl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(rgl,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_maf50.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(igl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(igl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(igl,mtra),lwd=3)


par(new=T)

data<-read.table("Pre_maf1.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(postc,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(postc,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(postc,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_mafFixed.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(dgl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(dgl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(dgl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,0.5,1),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Proportion", cex.lab = 1.5, line = 2)
put.fig.letter(label="c) proportions", location="topleft",cex=1.5)
#legend(17,1,legend=c(expression("AF "<=" 0.05"),expression("0.05 < AF "<=" 0.1"),expression("0.1 < AF "<=" 0.25"),
#						expression("0.25 < AF "<=" 0.5"),expression("0.5 < AF "<=" 1"),"AF = 1"),
#						col = transp(c(glc,fitn,rgl,igl,postc,dgl),.8),lwd=3,lty=1,box.lty=0,cex=1.0,pt.cex=2)

##############################################################################
##############################################################################
##############################################################################

data<-read.table("Pre_s01.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,100000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,100000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,100000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(glc,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,50000,100000),labels=c("0","50000","100000"),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Counts", cex.lab = 1.5, line = 2)
put.fig.letter(label="d) very weak mutations", location="topleft",cex=1.5)
legend(10,100000,legend=c(expression(plain("|")~italic('s')~plain("| ") <= " 0.001")),
						col = transp(c(glc),.8),lwd=3,lty=1,box.lty=0,cex=1.0,pt.cex=2)

#################################

data<-read.table("Pre_s10.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,40000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(rgl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,40000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(rgl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,40000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(rgl,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_s20.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,40000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(igl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,40000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(igl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,40000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(igl,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_s50.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,40000),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(dgl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,40000),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(dgl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,40000),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(dgl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,20000,40000),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Counts", cex.lab = 1.5, line = 2)
put.fig.letter(label="e) less weak mutations", location="topleft",cex=1.5)
legend(9,40000,legend=c(expression(plain('0.001 < ')~plain("|")~italic('s')~plain("| ") <= " 0.01"),
						expression(plain('0.01 < ')~plain("|")~italic('s')~plain("| ") <= " 0.1"),
						expression(plain("|")~italic('s')~plain("| > 0.1"))),
						col = transp(c(rgl,igl,dgl),.8),lwd=3,lty=1,box.lty=0,cex=1.0,pt.cex=2)

#############################################################

tots<-read.table("Total_prebott.txt")
data<-read.table("Pre_s01.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(sims,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(glc,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_s10.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(rgl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(rgl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(rgl,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_s20.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(igl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(igl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(igl,mtra),lwd=3)

par(new=T)

data<-read.table("Pre_s50.txt",h=F)
data<-data/tots
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',
	xlab="",main="",col=transp(dgl,tra), cex.lab = 1.5)
#rect(xleft=1,ybottom=-10,xright=6,ytop=30,col=transp(paste(highl[1]),rectr),border=NA)
for (i in 2:100){
par(new=T)
plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n', axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(dgl,tra))
}
par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n', axes=F,xlab="",main="",col=transp(dgl,mtra),lwd=3)
box(bty="l")
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(2,c(0,0.5,1),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Proportion", cex.lab = 1.5, line = 2)
put.fig.letter(label="f) proportions", location="topleft",cex=1.5)
#legend(17,1,legend=c(expression("AF "<=" 0.05"),expression("0.05 < AF "<=" 0.1"),expression("0.1 < AF "<=" 0.25"),
#						expression("0.25 < AF "<=" 0.5"),expression("0.5 < AF "<=" 1"),"AF = 1"),
#						col = transp(c(glc,fitn,rgl,igl,postc,dgl),.8),lwd=3,lty=1,box.lty=0,cex=1.0,pt.cex=2)

dev.off()
