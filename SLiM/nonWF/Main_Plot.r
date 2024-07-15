glc<-"#253741"
fitn<-"#DEC102" #"#FF7F00FF"
rgl<-"#FF7F00FF" #"#E41A1CFF"
dgl<-"#984EA3FF"
igl<-"#4DAF4AFF"
#postc<-"#E41A1CFF"
sims<-"grey60"
sims2<-"grey70"
sims3<-"grey50"
emp<-"tan1"
tra<-0.5
tra2<-0.075
rectr<-0.1
mtra<-0.8
fitn1<-"#99d8c9"
fitn2<-"#66c2a4"
fitn3<-"#41ae76"
fitn4<-"#238b45"
fitn5<-"#005824"
ext<-"#d95f02"
surv<-"#1b9e77"
pend<-"grey50"

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

pdf("nonWF_Final_New.pdf",width=9.5, height=6.5)
par(mfrow=c(2,3))

probs = c(1, 18/20, 3/20, 0, 0)
barplot(probs, space = 1, yaxt='n',axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(glc,tra))

box(bty="l")
axis(2,seq(0,1,by=0.2),labels=c(as.character(c(0,0.2,0.4,0.6,0.8,1))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(1,c(1.5,3.5,5.5,7.5,9.5),labels=c(as.character(c("50","100","250","500","1000"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
title(xlab = expression(italic('K')), cex.lab = 1.5, line = 2)
title(ylab = "Probability", cex.lab = 1.5, line = 2)
put.fig.letter(label="a) Extinction probabilities", location="topleft",cex=1.5)

# All input data here below are output from the various SLiM simulations.

data<-read.table("Fin_Bt100_sizes.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot((data$V1),type='l',ylim=c(0,350000),yaxt='n',axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
points(c(25-10.1,25-5.7,25-3.4,25-2.18,25-0),c(350,15000,120000,127000,225000),cex=2,pch=20,col=emp)
segments(c(25-10.1,25-5.7,25-3.4,25-2.18), c(350,15000,120000,127000), c(25-5.7,25-3.4,25-2.18,25-0), c(15000,120000,127000,225000), lwd=2, col=emp)
for (i in 2:100){
if (length(which(data[,i]%in%NA))!=0){
par(new=T)
plot(min((which(data[,i]%in%NA)))-2,-(data[min(which(data[,i]%in%NA))-1,i]),type='l',axes=F,col=transp(sims,tra),ylim=c(0,350000),main="",xlab="",ylab="",xaxt='n',yaxt='n')
} else {
par(new=T)
plot((data[,i]),type='l',ylim=c(0,350000),yaxt='n',axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}}
par(new=T)
plot(ms,type='l',ylim=c(0,350000),yaxt='n',ylab="",axes=F,xaxt='n',xlab="",main="",col=transp(glc,mtra),lwd=3)

box(bty="l")
axis(2,seq(0,300000,by=50000),labels=c(as.character(c(0,50,100,150,200,250,300))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .4, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = expression(italic('N')[c]*" (x 10"^3*")"), cex.lab = 1.5, line = 2)
put.fig.letter(label=expression(plain("b) Bottleneck ")~italic('K')~plain(" = 100")), location="topleft",cex=1.5)

#######

data<-read.table("Fin_Bt250_sizes.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot((data$V1),type='l',ylim=c(0,350000),yaxt='n',axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
points(c(25-10.1,25-5.7,25-3.4,25-2.18,25-0),c(350,15000,120000,127000,225000),cex=2,pch=20,col=emp)
segments(c(25-10.1,25-5.7,25-3.4,25-2.18), c(350,15000,120000,127000), c(25-5.7,25-3.4,25-2.18,25-0), c(15000,120000,127000,225000), lwd=2, col=emp)
for (i in 2:100){
if (length(which(data[,i]%in%NA))!=0){
par(new=T)
plot(min((which(data[,i]%in%NA)))-2,-(data[min(which(data[,i]%in%NA))-1,i]),type='l',axes=F,col=transp(sims,tra),ylim=c(0,350000),main="",xlab="",ylab="",xaxt='n',yaxt='n')
} else {
par(new=T)
plot((data[,i]),type='l',ylim=c(0,350000),yaxt='n',axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}}
par(new=T)
plot(ms,type='l',ylim=c(0,350000),yaxt='n',ylab="",axes=F,xaxt='n',xlab="",main="",col=transp(glc,mtra),lwd=3)

box(bty="l")
axis(2,seq(0,300000,by=50000),labels=c(as.character(c(0,50,100,150,200,250,300))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .4, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = expression(italic('N')[c]*" (x 10"^3*")"), cex.lab = 1.5, line = 2)
put.fig.letter(label=expression(plain("c) Bottleneck ")~italic('K')~plain(" = 250")), location="topleft",cex=1.5)

#######

data<-read.table("Fin_Bt500_sizes.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot((data$V1),type='l',ylim=c(0,350000),yaxt='n',axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
points(c(25-10.1,25-5.7,25-3.4,25-2.18,25-0),c(350,15000,120000,127000,225000),cex=2,pch=20,col=emp)
segments(c(25-10.1,25-5.7,25-3.4,25-2.18), c(350,15000,120000,127000), c(25-5.7,25-3.4,25-2.18,25-0), c(15000,120000,127000,225000), lwd=2, col=emp)
for (i in 2:100){
if (length(which(data[,i]%in%NA))!=0){
par(new=T)
plot(min((which(data[,i]%in%NA)))-2,-(data[min(which(data[,i]%in%NA))-1,i]),type='l',axes=F,col=transp(sims,tra),ylim=c(0,350000),main="",xlab="",ylab="",xaxt='n',yaxt='n')
} else {
par(new=T)
plot((data[,i]),type='l',ylim=c(0,350000),yaxt='n',axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}}
par(new=T)
plot(ms,type='l',ylim=c(0,350000),yaxt='n',ylab="",axes=F,xaxt='n',xlab="",main="",col=transp(glc,mtra),lwd=3)

box(bty="l")
axis(2,seq(0,300000,by=50000),labels=c(as.character(c(0,50,100,150,200,250,300))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .4, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = expression(italic('N')[c]*" (x 10"^3*")"), cex.lab = 1.5, line = 2)
put.fig.letter(label=expression(plain("d) Bottleneck ")~italic('K')~plain(" = 500")), location="topleft",cex=1.5)

#######

data<-read.table("Fin_Bt1000_sizes.txt",h=F)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
plot((data$V1),type='l',ylim=c(0,350000),yaxt='n',axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
points(c(25-10.1,25-5.7,25-3.4,25-2.18,25-0),c(350,15000,120000,127000,225000),cex=2,pch=20,col=emp)
segments(c(25-10.1,25-5.7,25-3.4,25-2.18), c(350,15000,120000,127000), c(25-5.7,25-3.4,25-2.18,25-0), c(15000,120000,127000,225000), lwd=2, col=emp)
for (i in 2:100){
if (length(which(data[,i]%in%NA))!=0){
par(new=T)
plot(min((which(data[,i]%in%NA)))-2,-(data[min(which(data[,i]%in%NA))-1,i]),type='l',axes=F,col=transp(sims,tra),ylim=c(0,350000),main="",xlab="",ylab="",xaxt='n',yaxt='n')
} else {
par(new=T)
plot((data[,i]),type='l',ylim=c(0,350000),yaxt='n',axes=F,ylab="",xaxt='n',xlab="",main="",col=transp(sims,tra))
}}
par(new=T)
plot(ms,type='l',ylim=c(0,350000),yaxt='n',ylab="",axes=F,xaxt='n',xlab="",main="",col=transp(glc,mtra),lwd=3)

box(bty="l")
axis(2,seq(0,300000,by=50000),labels=c(as.character(c(0,50,100,150,200,250,300))),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .4, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = expression(italic('N')[c]*" (x 10"^3*")"), cex.lab = 1.5, line = 2)
put.fig.letter(label=expression(plain("e) Bottleneck ")~italic('K')~plain(" = 1000")), location="topleft",cex=1.5)

data<-read.table("Fin_Bt50_RGL.txt",h=F)
data<-exp(data)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
#plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#for (i in 2:100){
#if (length(which(data[,i]%in%NA))!=0){
#par(new=T)
#plot(min((which(data[,i]%in%NA)))-2,-(data[min(which(data[,i]%in%NA))-1,i]),axes=F,type='l',col=transp(sims,tra),xlim=c(0,24),ylim=c(0,8),main="",xlab="",ylab="",xaxt='n',yaxt='n')
#} else {
#par(new=T)
#plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#}}
#par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',xlab="",axes=F,main="",col=transp(fitn1,mtra),lwd=3)
par(new=T)

data<-read.table("Fin_Bt100_RGL.txt",h=F)
data<-exp(data)
ms<-apply(data[,-16],1,function(x) mean(x,na.rm=T))
#plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#for (i in 2:100){
#if (length(which(data[,i]%in%NA))!=0){
#par(new=T)
#plot(min((which(data[,i]%in%NA)))-2,-(data[min(which(data[,i]%in%NA))-1,i]),axes=F,type='l',col=transp(sims,tra),xlim=c(0,24),ylim=c(0,8),main="",xlab="",ylab="",xaxt='n',yaxt='n')
#} else {
#par(new=T)
#plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#}}
#par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',xlab="",axes=F,main="",col=transp(fitn2,mtra),lwd=3)
par(new=T)

data<-read.table("Fin_Bt250_RGL.txt",h=F)
data<-exp(data)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
#plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#for (i in 2:100){
#if (length(which(data[,i]%in%NA))!=0){
#par(new=T)
#plot(min((which(data[,i]%in%NA)))-2,-(data[min(which(data[,i]%in%NA))-1,i]),axes=F,type='l',col=transp(sims,tra),xlim=c(0,24),ylim=c(0,8),main="",xlab="",ylab="",xaxt='n',yaxt='n')
#} else {
#par(new=T)
#plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#}}
#par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',xlab="",axes=F,main="",col=transp(fitn3,mtra),lwd=3)
par(new=T)

data<-read.table("Fin_Bt500_RGL.txt",h=F)
data<-exp(data)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
#plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#for (i in 2:100){
#if (length(which(data[,i]%in%NA))!=0){
#par(new=T)
#plot(min((which(data[,i]%in%NA)))-2,-(data[min(which(data[,i]%in%NA))-1,i]),axes=F,type='l',col=transp(sims,tra),xlim=c(0,24),ylim=c(0,8),main="",xlab="",ylab="",xaxt='n',yaxt='n')
#} else {
#par(new=T)
#plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#}}
#par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',xlab="",axes=F,main="",col=transp(fitn4,mtra),lwd=3)
par(new=T)

data<-read.table("Fin_Bt1000_RGL.txt",h=F)
data<-exp(data)
ms<-apply(data,1,function(x) mean(x,na.rm=T))
#plot(0:24,(data$V1),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#for (i in 2:100){
#if (length(which(data[,i]%in%NA))!=0){
#par(new=T)
#plot(min((which(data[,i]%in%NA)))-2,-(data[min(which(data[,i]%in%NA))-1,i]),axes=F,type='l',col=transp(sims,tra),xlim=c(0,24),ylim=c(0,8),main="",xlab="",ylab="",xaxt='n',yaxt='n')
#} else {
#par(new=T)
#plot(0:24,(data[,i]),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',axes=F,xlab="",main="",col=transp(sims,tra))
#}}
#par(new=T)
plot((ms),type='l',ylim=c(0,1),yaxt='n',ylab="",xaxt='n',xlab="",axes=F,main="",col=transp(fitn5,mtra),lwd=3)
box(bty="l")
axis(2,seq(0,1,by=0.2),cex.axis=1.5,tck=-0.005,mgp=c(3, .3, 0))
axis(1,c(2,7,12,17,22),labels=c(as.character(c("0","5","10","15","20"))),cex.axis=1.5,tck=-0.005,mgp=c(3, .4, 0))
abline(v=2,lty=5,col="#E41A1CFF")
title(xlab = "Generations", cex.lab = 1.5, line = 2)
title(ylab = "Fitness", cex.lab = 1.5, line = 2)
points(c(25-10.1,25-5.7,25-3.4,25-2.18,25-0),c(350,15000,120000,127000,225000),cex=2,pch=20,col='red')
segments(c(25-10.1,25-5.7,25-3.4,25-2.18), c(350,15000,120000,127000), c(25-5.7,25-3.4,25-2.18,25-0), c(15000,120000,127000,225000), lwd=2, col='red')
put.fig.letter(label="f) Fitness", location="topleft",cex=1.5)

legend(18,1.05,legend=c("50","100","250","500","1000"),col = transp(c(fitn1,fitn2,fitn3,fitn4,fitn5),mtra),
		lwd=3,lty=1,box.lty=0,cex=1.0,pt.cex=2)

#################

dev.off()













