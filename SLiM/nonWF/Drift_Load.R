#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

files<-read.table("file_list.txt",h=F)
if (dim(files)[1]==25){

Nes<-read.table(args[1],h=F)
DL<-c() # drift load

data<-read.table(paste(files[1,1]),h=F)
Ne<-Nes[1,1]
colnames(data)<-c("id1","ID","mut","pos","s","h","pop","GEN","N")

drift<-data[data$N==(2*Ne),]
ids<-drift$ID
DL[1]<-sum(drift$s)

for (i in 2:length(files[,1])){

data<-read.table(paste(files[i,1]),h=F)
Ne<-Nes[i,1]
colnames(data)<-c("id1","ID","mut","pos","s","h","pop","GEN","N")

drift<-data[data$N==(2*Ne),]
#ind<-which(drift$ID%in%ids==FALSE)
#ids<-sort(c(ids,drift$ID[ind]))
#drift<-drift[ind,]
DL[i]<-sum(drift$s)
}

write.table(DL,"DL.txt",quote=F,col.names=F,row.names=F,sep="\t")

} else {

Nes<-read.table(args[1],h=F)
DL<-c() # drift load

data<-read.table(paste(files[1,1]),h=F)
Ne<-Nes[1,1]
colnames(data)<-c("id1","ID","mut","pos","s","h","pop","GEN","N")

drift<-data[data$N==(2*Ne),]
ids<-drift$ID
DL[1]<-sum(drift$s)

for (i in 1:length(files[,1])){

data<-read.table(paste(files[i,1]),h=F)
Ne<-Nes[i,1]
colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")

drift<-data[data$N==(2*Ne),]
#ind<-which(drift$ID%in%ids==FALSE)
#ids<-sort(c(ids,drift$ID[ind]))
#drift<-drift[ind,]
DL[i]<-sum(drift$s)
}

DL[(length(files[,1])+1):25]<-NA

write.table(DL,"DL.txt",quote=F,col.names=F,row.names=F,sep="\t")
}
