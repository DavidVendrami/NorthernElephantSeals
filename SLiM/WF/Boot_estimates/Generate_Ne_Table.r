# Generate Ne tables

data <- read.table("Inputs.txt",h=F)
cnam<-read.table("Input_colnames.txt",h=F)
colnames(data)<-t(cnam[,1])

out <- matrix( ,ncol=100,nrow=25)

for (i in 1:100){
out[1,i]<-round(data$LGM[i] * exp(data$GRLGM[i]*(data$YYYYY[i]-data$XXXXX[i])))
out[(2:7),i]<-data$BOTT[i]

for (j in 0:17){
out[(j+8),i]<-round(data$BOTT[i] * exp(data$GRBOT[i]*j))
}
}
write.table(out,"RAD_Boot_Ne.txt",quote=F,col.names=F,row.names=F,sep="\t")
