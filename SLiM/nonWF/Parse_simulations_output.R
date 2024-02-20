#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

files<-read.table("file_list.txt",h=F)
if (dim(files)[1]==25){

Nes<-read.table(args[1],h=F)
GL<-c() # Genetic load
RGL<-c() # Realized genetic load
MGL<-c() # Masked genetic load
LD<-c() # Number of mildly deleterious mutations
HD<-c() # Number of highly deleterious mutations
fLD<-c() # Number of fixed mildly deleterious mutations
fHD<-c() # Number of fixed highly deleterious mutations


for (i in 1:length(files[,1])){

data<-read.table(paste(files[i,1]),h=F)
Ne<-Nes[i,1]
colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")

af<-data$N/(2*Ne)
data$q<-af
gl<-data$s*data$q
data$gl<-gl
GL[i]<-sum(data$gl)


rgl_p1<-(data$q^2)*data$s
rgl_p2<-data$q*(1-data$q)*data$h*data$s
srgl_p1<-sum(rgl_p1)
srgl_p2<-2*sum(rgl_p2)
RGL[i]<-srgl_p1+srgl_p2

MGL[i]<-GL[i]-RGL[i]

LD[i]<-length(which(data$mut=="m2"))

HD[i]<-length(which(data$mut=="m3"))

fLD[i]<-length(which(data$mut=="m2"&data$N==(Ne*2)))

fHD[i]<-length(which(data$mut=="m3"&data$N==(Ne*2)))
}

write.table(GL,"GL.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(RGL,"RGL.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(MGL,"MGL.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(LD,"LD.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(HD,"HD.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(fLD,"fLD.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(fHD,"fHD.txt",quote=F,col.names=F,row.names=F,sep="\t")

} else {

Nes<-read.table(args[1],h=F)
GL<-c() # Genetic load
RGL<-c() # Realized genetic load
MGL<-c() # Masked genetic load
LD<-c() # Number of mildly deleterious mutations
HD<-c() # Number of highly deleterious mutations
fLD<-c() # Number of fixed mildly deleterious mutations
fHD<-c() # Number of fixed highly deleterious mutations


for (i in 1:length(files[,1])){

data<-read.table(paste(files[i,1]),h=F)
Ne<-Nes[i,1]
colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")

af<-data$N/(2*Ne)
data$q<-af
gl<-data$s*data$q
data$gl<-gl
GL[i]<-sum(data$gl)


rgl_p1<-(data$q^2)*data$s
rgl_p2<-data$q*(1-data$q)*data$h*data$s
srgl_p1<-sum(rgl_p1)
srgl_p2<-2*sum(rgl_p2)
RGL[i]<-srgl_p1+srgl_p2

MGL[i]<-GL[i]-RGL[i]

LD[i]<-length(which(data$mut=="m2"))

HD[i]<-length(which(data$mut=="m3"))

fLD[i]<-length(which(data$mut=="m2"&data$N==(Ne*2)))

fHD[i]<-length(which(data$mut=="m3"&data$N==(Ne*2)))
}

GL[(length(files[,1])+1):25]<-NA
RGL[(length(files[,1])+1):25]<-NA
MGL[(length(files[,1])+1):25]<-NA
LD[(length(files[,1])+1):25]<-NA
HD[(length(files[,1])+1):25]<-NA
fLD[(length(files[,1])+1):25]<-NA
fHD[(length(files[,1])+1):25]<-NA

write.table(GL,"GL.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(RGL,"RGL.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(MGL,"MGL.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(LD,"LD.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(HD,"HD.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(fLD,"fLD.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(fHD,"fHD.txt",quote=F,col.names=F,row.names=F,sep="\t")}
