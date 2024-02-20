## Calculate sMLH and g2 for msats and SNPs.
library(dplyr)
library(inbreedR)
library(tidyr)

################ Microsatellites
msats<-read.table("msats_genotypes.txt",h=F,row.names=1)
msats[msats==0]<-NA
outm<-read.table("msats_dataset_description.txt",h=T)
outm$Animal_ID <- gsub("-","",outm$Animal_ID)
msats<-convert_raw(msats)

m_smlh<-sMLH(msats)
m_smlh<-data.frame(Animal_ID=row.names(msats),sMLH=m_smlh)
out_m<-left_join(outm,m_smlh,by="Animal_ID")

write.table(out_m,"msats_Samples_sMLH",quote=F,col.names=T,row.names=F,sep="\t")

g2_msats_100 <- g2_microsats(msats, nperm = 100, nboot = 100, CI = 0.95)
g2_msats_1000 <- g2_microsats(msats, nperm = 100, nboot = 1000, CI = 0.95)
g2_msats_10000 <- g2_microsats(msats, nperm = 100, nboot = 10000, CI = 0.95)

write.table(g2_msats_100$g2_boot,"g2/g2msats_100.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(g2_msats_1000$g2_boot,"g2/g2msats_1000.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(g2_msats_10000$g2_boot,"g2/g2msats_10000.txt",quote=F,col.names=F,row.names=F,sep="\t")

################ SNPs
raw<-read.table("NES_notImp.raw",h=T) # Genotype data
out<-read.table("74_Samples_sMLH.txt",h=T) # This is the file where sMLH info will be saved (currently an empty file)

# Rename 'raw' IDs to match those of 'out'
iid<-raw$IID
iid<-gsub('-','',iid)
iid<-gsub('_.*','',iid)
raw$IID<-iid

# Calculate genome-wide sMLH
geno<-raw[,-c(1:6)] # Remove columns with no genotypic data
geno[geno==2]<-0
smlh<-sMLH(geno)

# Put it onto the table
temp1<-data.frame(ID=raw$IID,sMLH=smlh)
temp2<-data.frame(ID=out$ID)
temp3<-left_join(temp2,temp1,by="ID")

out$sMLH_GenomeWide <- temp3$sMLH

# sMLH separately for each chromosome
cols<-colnames(geno)
cols<-gsub('.[0-9]*_[ACTG]','',cols) # Keep only chr name
cols[12778:15051]<-"Unplaced"
fac<-factor(cols, levels = c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5",
	"HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12",
	"HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","Unplaced"))
index<-table(fac) 

for (i in 1:length(index)){
if (i==1){
temp<-geno[,1:index[i]]
assign(paste("chr_",i,sep=""),sMLH(temp))
} else {
temp<-geno[,(((index[i-1]+sum(index[1:(i-2)]))+1):(index[i]+sum(index[1:(i-1)])))]
assign(paste("chr_",i,sep=""),sMLH(temp))
}
}

out$sMLH_Chr1 <- chr_1
out$sMLH_Chr2 <- chr_2
out$sMLH_Chr3 <- chr_3
out$sMLH_Chr4 <- chr_4
out$sMLH_Chr5 <- chr_5
out$sMLH_Chr6 <- chr_6
out$sMLH_Chr7 <- chr_7
out$sMLH_Chr8 <- chr_8
out$sMLH_Chr9 <- chr_9
out$sMLH_Chr10 <- chr_10
out$sMLH_Chr11 <- chr_11
out$sMLH_Chr12 <- chr_12
out$sMLH_Chr13 <- chr_13
out$sMLH_Chr14 <- chr_14
out$sMLH_Chr15 <- chr_15
out$sMLH_Chr16 <- chr_16
out$sMLH_Chr17 <- chr_17

# Now let's generate sliding window sMLH alogn MHC chromosome (i.e.: HiC_scaffold_16, length: 135,623,083).
# Chr 16 has 675 SNPs.
# Let's use a window size of 1Mb and a step of 100kb (based on LD decay plot).
stepw <- 100000
size <- 1000000
# Let's first create a "map" file
cols<-colnames(raw)[7:dim(raw)[2]]
cols<-gsub("_[ACTG]","",cols)
map<-data.frame(cols=cols)
map<-separate_wider_delim(map, cols = cols, delim = ".", names = c("CHR", "POS"))
map$POS<-as.numeric(map$POS)

# Get starting positions of windows
starts <- seq(1,135623083,by=stepw)



slw_report<-matrix( ,nrow=length(starts), ncol=4)
slw_smlh<-matrix( ,nrow=74, ncol=length(starts))
for (i in 1:length(starts)){
	if (starts[i]+size <= 135623083) {
		ind<-which(map$CHR=='HiC_scaffold_16' & map$POS >= starts[i] & map$POS <= (starts[i]+size))
		slw_report[i,1]<-"HiC_scaffold_16"
		slw_report[i,2]<-starts[i]
		slw_report[i,3]<-starts[i]+size
		slw_report[i,4]<-length(ind)
		tempg<-geno[,ind]
		slw_smlh[,i]<-sMLH(tempg)
	} else {
		ind<-which(map$CHR=='HiC_scaffold_16' & map$POS >= starts[i] & map$POS <= 135623083)
		slw_report[i,1]<-"HiC_scaffold_16"
		slw_report[i,2]<-starts[i]
		slw_report[i,3]<-135623083
		slw_report[i,4]<-length(ind)
		tempg<-geno[,ind]
		slw_smlh[,i]<-sMLH(tempg)
}
}


ID=seq(1,1357,by=1)
sl<-c()
sl[1:1357]<-"Win_"
win_ID<-paste(sl,ID,sep="")
sl_rep_out<-data.frame(Window_ID=win_ID,CHR=slw_report[,1],Start=slw_report[,2],Stop=slw_report[,3],N_SNPs=slw_report[,4])
sl_smlh_out<-data.frame(Sample_ID=raw$IID,slw_smlh)

write.table(sl_rep_out,"Sliding_Windows_Description.txt",quote=F,col.names=T,row.names=F,sep="\t")
write.table(sl_smlh_out,"Sliding_Windows_sMLH.txt",quote=F,col.names=T,row.names=F,sep="\t")


# g2 calculations, with different amounts of bootstrap from plotting purposes
g2_100<-g2_snps(geno,nperm=100,nboot=100,CI=0.95,parallel=F,ncores=NULL)
g2_1000<-g2_snps(geno,nperm=100,nboot=1000,CI=0.95,parallel=F,ncores=NULL)
g2_10000<-g2_snps(geno,nperm=100,nboot=10000,CI=0.95,parallel=F,ncores=NULL)

write.table(g2_100$g2_boot,"g2/g2snps_100.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(g2_1000$g2_boot,"g2/g2snps_1000.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(g2_10000$g2_boot,"g2/g2snps_10000.txt",quote=F,col.names=F,row.names=F,sep="\t")
