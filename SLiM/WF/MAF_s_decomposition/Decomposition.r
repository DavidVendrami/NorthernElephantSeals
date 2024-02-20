### Purging ###

# Setwd within one of the rep folders containing SLiM output
out_totpre<-matrix( ,ncol=100,nrow=25)
out_totnew<-matrix( ,ncol=100,nrow=25)
out_pre_maf05<-matrix( ,ncol=100,nrow=25)
out_pre_maf10<-matrix( ,ncol=100,nrow=25)
out_pre_maf25<-matrix( ,ncol=100,nrow=25)
out_pre_maf50<-matrix( ,ncol=100,nrow=25)
out_pre_maf1<-matrix( ,ncol=100,nrow=25)
out_pre_fixed<-matrix( ,ncol=100,nrow=25)
out_pre_s01<-matrix( ,ncol=100,nrow=25)
out_pre_s10<-matrix( ,ncol=100,nrow=25)
out_pre_s20<-matrix( ,ncol=100,nrow=25)
out_pre_s50<-matrix( ,ncol=100,nrow=25)
out_new_maf05<-matrix( ,ncol=100,nrow=25)
out_new_maf10<-matrix( ,ncol=100,nrow=25)
out_new_maf25<-matrix( ,ncol=100,nrow=25)
out_new_maf50<-matrix( ,ncol=100,nrow=25)
out_new_maf1<-matrix( ,ncol=100,nrow=25)
out_new_fixed<-matrix( ,ncol=100,nrow=25)
out_new_s01<-matrix( ,ncol=100,nrow=25)
out_new_s10<-matrix( ,ncol=100,nrow=25)
out_new_s20<-matrix( ,ncol=100,nrow=25)
out_new_s50<-matrix( ,ncol=100,nrow=25)
out_pre_gl05<-matrix( ,ncol=100,nrow=25)
out_new_gl05<-matrix( ,ncol=100,nrow=25)
out_pre_gl10<-matrix( ,ncol=100,nrow=25)
out_new_gl10<-matrix( ,ncol=100,nrow=25)
out_pre_gl25<-matrix( ,ncol=100,nrow=25)
out_new_gl25<-matrix( ,ncol=100,nrow=25)
out_pre_gl50<-matrix( ,ncol=100,nrow=25)
out_new_gl50<-matrix( ,ncol=100,nrow=25)

for (j in 1:100){
setwd(paste("../rep_",j,sep=""))

files<-read.table("file_list.txt",h=F)

tot_pre<-c()
tot_new<-c()
pre_maf05<-c()
pre_maf10<-c()
pre_maf25<-c()
pre_maf50<-c()
new_maf05<-c()
new_maf10<-c()
new_maf25<-c()
new_maf50<-c()
pre_maf1<-c()
new_maf1<-c()
pre_fixed<-c()
new_fixed<-c()
pre_s01<-c()
pre_s10<-c()
pre_s20<-c()
pre_s50<-c()
new_s01<-c()
new_s10<-c()
new_s20<-c()
new_s50<-c()
pre_gl05<-c()
new_gl05<-c()
pre_gl10<-c()
new_gl10<-c()
pre_gl25<-c()
new_gl25<-c()
pre_gl50<-c()
new_gl50<-c()

# tots
preb<-read.table(paste(files[1,1]),h=F)
colnames(preb)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
pre_mut<-preb$id2

for (i in 1:length(files[,1])){
data<-read.table(paste(files[i,1]),h=F)
colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
tot_pre<-c(tot_pre,length(which(data$id2%in%pre_mut)))
tot_new<-c(tot_new,length(which(!data$id2%in%pre_mut)))
}

# mafs
Nes<-read.table("/prj/furseal-genome/David/NES/SLiM/SLiM_3.7.1/WF/Ne_RAD_table.txt",h=T)
preb<-read.table(paste(files[1,1]),h=F)
Ne<-Nes[1,2]
colnames(preb)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
pre_mut<-preb$id2

for (i in 1:length(files[,1])){
data<-read.table(paste(files[i,1]),h=F)
Ne<-Nes[i,2]
colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
data$q<-data$N/(2*Ne)
pre_maf05<-c(pre_maf05,length(which(data$id2[data$q<=0.05] %in% pre_mut)))
new_maf05<-c(new_maf05,length(which(!(data$id2[data$q<=0.05] %in% pre_mut))))

pre_maf10<-c(pre_maf10,length(which(data$id2[data$q>0.05 & data$q<=0.1] %in% pre_mut)))
new_maf10<-c(new_maf10,length(which(!(data$id2[data$q>0.05 & data$q<=0.1] %in% pre_mut))))

pre_maf25<-c(pre_maf25,length(which(data$id2[data$q>0.1 & data$q<=0.25] %in% pre_mut)))
new_maf25<-c(new_maf25,length(which(!(data$id2[data$q>0.1 & data$q<=0.25] %in% pre_mut))))

pre_maf50<-c(pre_maf50,length(which(data$id2[data$q>0.25 & data$q<=0.5] %in% pre_mut)))
new_maf50<-c(new_maf50,length(which(!(data$id2[data$q>0.25 & data$q<=0.5] %in% pre_mut))))

pre_maf1<-c(pre_maf1,length(which(data$id2[data$q>0.5 & data$q<1] %in% pre_mut)))
new_maf1<-c(new_maf1,length(which(!(data$id2[data$q>0.5 & data$q<1] %in% pre_mut))))

pre_fixed<-c(pre_fixed,length(which(data$id2[data$q==1] %in% pre_mut)))
new_fixed<-c(new_fixed,length(which(!(data$id2[data$q==1] %in% pre_mut))))
}

# s
preb<-read.table(paste(files[1,1]),h=F)
colnames(preb)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
pre_mut<-preb$id2

for (i in 1:length(files[,1])){
data<-read.table(paste(files[i,1]),h=F)
Ne<-Nes[i,2]
colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
data$q<-data$N/(2*Ne)
data$s<-data$s*(-1)
pre_s01<-c(pre_s01,length(which(data$id2[data$s<=0.001] %in% pre_mut)))
new_s01<-c(new_s01,length(which(!(data$id2[data$s<=0.001] %in% pre_mut))))

pre_s10<-c(pre_s10,length(which(data$id2[data$s>0.001 & data$s<=0.01] %in% pre_mut)))
new_s10<-c(new_s10,length(which(!(data$id2[data$s>0.001 & data$s<=0.01] %in% pre_mut))))

pre_s20<-c(pre_s20,length(which(data$id2[data$s>0.01 & data$s<=0.1] %in% pre_mut)))
new_s20<-c(new_s20,length(which(!(data$id2[data$s>0.01 & data$s<=0.1] %in% pre_mut))))

pre_s50<-c(pre_s50,length(which(data$id2[data$s>0.1] %in% pre_mut)))
new_s50<-c(new_s50,length(which(!(data$id2[data$s>0.1] %in% pre_mut))))
}

# GLs
Nes<-read.table("/prj/furseal-genome/David/NES/SLiM/SLiM_3.7.1/WF/Ne_RAD_table.txt",h=T)
preb<-read.table(paste(files[1,1]),h=F)
colnames(preb)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
pre_mut<-preb$id2

for (i in 1:length(files[,1])){
data<-read.table(paste(files[i,1]),h=F)
Ne<-Nes[i,2]
colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
data$q<-data$N/(2*Ne)
data$gl<-(-data$s*data$q)
pre_gl05<-c(pre_gl05,length(which(data$id2[data$gl<=0.0001] %in% pre_mut)))
new_gl05<-c(new_gl05,length(which(!(data$id2[data$gl<=0.0001] %in% pre_mut))))

pre_gl10<-c(pre_gl10,length(which(data$id2[data$gl>0.0001 & data$gl<=0.001] %in% pre_mut)))
new_gl10<-c(new_gl10,length(which(!(data$id2[data$gl>0.0001 & data$gl<=0.001] %in% pre_mut))))

pre_gl25<-c(pre_gl25,length(which(data$id2[data$gl>0.001 & data$gl<=0.01] %in% pre_mut)))
new_gl25<-c(new_gl25,length(which(!(data$id2[data$gl>0.001 & data$gl<=0.01] %in% pre_mut))))

pre_gl50<-c(pre_gl50,length(which(data$id2[data$gl>0.01] %in% pre_mut)))
new_gl50<-c(new_gl50,length(which(!(data$id2[data$gl>0.01] %in% pre_mut))))

}





out_totpre[,j]<-tot_pre
out_totnew[,j]<-tot_new
out_pre_maf05[,j]<-pre_maf05
out_pre_maf10[,j]<-pre_maf10
out_pre_maf25[,j]<-pre_maf25
out_pre_maf50[,j]<-pre_maf50
out_pre_maf1[,j]<-pre_maf1
out_pre_fixed[,j]<-pre_fixed
out_pre_s01[,j]<-pre_s01
out_pre_s10[,j]<-pre_s10
out_pre_s20[,j]<-pre_s20
out_pre_s50[,j]<-pre_s50
out_new_maf05[,j]<-new_maf05
out_new_maf10[,j]<-new_maf10
out_new_maf25[,j]<-new_maf25
out_new_maf50[,j]<-new_maf50
out_new_maf1[,j]<-new_maf1
out_new_fixed[,j]<-new_fixed
out_new_s01[,j]<-new_s01
out_new_s10[,j]<-new_s10
out_new_s20[,j]<-new_s20
out_new_s50[,j]<-new_s50
out_pre_gl05[,j]<-pre_gl05
out_new_gl05[,j]<-new_gl05
out_pre_gl10[,j]<-pre_gl10
out_new_gl10[,j]<-new_gl10
out_pre_gl25[,j]<-pre_gl25
out_new_gl25[,j]<-new_gl25
out_pre_gl50[,j]<-pre_gl50
out_new_gl50[,j]<-new_gl50
}

write.table(out_totpre,"../Total_prebott.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_totnew,"../Total_new.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_maf05,"../Pre_maf05.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_maf10,"../Pre_maf10.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_maf25,"../Pre_maf25.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_maf50,"../Pre_maf50.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_maf1,"../Pre_maf1.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_fixed,"../Pre_mafFixed.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_s01,"../Pre_s01.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_s10,"../Pre_s10.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_s20,"../Pre_s20.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_s50,"../Pre_s50.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_maf05,"../New_maf05.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_maf10,"../New_maf10.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_maf25,"../New_maf25.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_maf50,"../New_maf50.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_maf1,"../New_maf1.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_fixed,"../New_mafFixed.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_s01,"../New_s01.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_s10,"../New_s10.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_s20,"../New_s20.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_s50,"../New_s50.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_gl05,"../PreGLmin.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_gl10,"../PreGLweak.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_gl25,"../PreGLmod.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_pre_gl50,"../PreGLstr.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_gl05,"../NewGLmin.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_gl10,"../NewGLweak.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_gl25,"../NewGLmod.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(out_new_gl50,"../NewGLstr.txt",quote=F,col.names=F,row.names=F,sep="\t")
