################################################################
### Let's calcualte heterozygosity from genotype likelihoods ###
################################################################

## Total Heterozygosity
# Seems to be  abit more laborious than expected... It requires to generate single samples SFS and use that to calculate heterozygosity
# Samples' bam files are listed in 'bam_files.txt'
# Let's then loop through it to generate single samples SFS
cat bam_files.txt | while read i
do
angsd -i $i -anc Mirounga_angustirostris_HiC.fasta -ref Mirounga_angustirostris_HiC.fasta -uniqueOnly 1 -remove_bads 1 -baq 1 -minMapQ 20 -minQ 30 -doSaf 1 -P 8 -GL 1 -out ./ANGSD_het/${i%_CoordSorted_RmDup_NoSex.bam}
realSFS ./ANGSD_het/${i%_CoordSorted_RmDup_NoSex.bam}.saf.idx -maxIter 100 -fold 1 -P 8 > ./ANGSD_het/${i%_CoordSorted_RmDup_NoSex.bam}.ml
done

cd ANGSD_het/
ls -1 *.ml > samples_ml_list.txt
Rscript Calc_het.r

# Where 'calc_het.r' is:
# sams <- read.table("samples_ml_list.txt", h = F)
# 
# hets <- c()
# ids <- c()
# 
# for (i in sams[,1]){
# ids <- c(ids, paste(i))
# a <- scan(paste(i))
# hets <- c(hets, a[2]/sum(a))
# }
# 
# out <- data.frame(Sample_ID = ids, GL_Heterozygosity = hets)
# write.table(out, "GL_Heterozugosities.txt", quote = F, col.names = T, row.names = F, sep = "\t")

########################################

## Let's now process Heterozygosity output:
library(dplyr)

data<-read.table("GL_Heterozugosities.txt",h=T) # Empty matrix except sample IDs in forst column to store GL based heterozygsity values.
# data<-data[1:96,]
data$Sample_ID<-gsub(".ml","",data$Sample_ID)
data$Sample_ID<-gsub("sample_","",data$Sample_ID)

conv<-read.table("IDs_conversion.txt",h=T)
colnames(data)[1]<-"barcode"
colnames(conv)[2]<-"Sample_ID"

temp<-left_join(data,conv,by="barcode")
ccs<-read.table("Phen.txt",h=T) # File containing phenotype information
ccs$Animal_ID <- gsub("ES","ES-",ccs$Animal_ID)
colnames(ccs)[1]<-"Sample_ID"
temp2<-left_join(temp,ccs,by="Sample_ID") 
out<-data.frame(Sample_ID=temp2$Sample_ID,Barcode=temp2$barcode,Class=temp2$Primary_classification, GL_Het=temp2$GL_Heterozygosity) # 1 = Worm
write.table(out,"GL_Heterozygosity_dataset.txt",quote=F,col.names=T,row.names=F,sep="\t")

##############################################
## Per chromosome heterozygosity
cat CHRs.txt | while read CHR # One chromosome ID per line
do
cat bam_files.txt | while read i
do
realSFS ./${i%_CoordSorted_RmDup.bam}.saf.idx -maxIter 100 -fold 1 -P 8 -r $CHR > ANGSD_GLhet_$CHR/${i%_CoordSorted_RmDup.bam}.ml
done

cd ANGSD_GLhet_$CHR
ls -1 *.ml > samples_ml_list.txt
Rscript Calc_het.r
rm *.ml
cd ../

done

for i in {1..17}
do 
mv ANGSD_GLhet_HiC_scaffold_$i/GL_Heterozugosities.txt ANGSD_GLhet_HiC_scaffolds/GL_Heterozugosities_$i.txt
done

cd ANGSD_GLhet_HiC_scaffolds
ls -1 GL_Heterozugosities_*.txt > GL_CHR_files_list.txt

########################################

## Let's now process Heterozygosity output:
library(dplyr)

conv<-read.table("IDs_conversion.txt",h=T)
out<-read.table('GL_Heterozygosity_dataset.txt',h=T) # Copy of GL_Heterozygosity_dataset.txt with 17 extra empty columns, one per chr.
ccs<-read.table("Phen.txt",h=T)
ccs$Animal_ID <- gsub("ES","ES-",ccs$Animal_ID)
colnames(ccs)[1]<-"Sample_ID"
colnames(conv)[2]<-"Sample_ID"

for (i in 1:17){

data<-read.table(paste('GL_Heterozugosities_',i,'.txt',sep=""),h=T)
# data<-data[1:96,]
data$Sample_ID<-gsub(".ml","",data$Sample_ID)
data$Sample_ID<-gsub("sample_","",data$Sample_ID)

colnames(data)[1]<-"barcode"

temp<-left_join(data,conv,by="barcode")
temp2<-left_join(temp,ccs,by="Sample_ID") 

assign(paste('chr_',i,sep=""),temp2$GL_Heterozygosity)
}

out$GL_het_Chr1 <- chr_1
out$GL_het_Chr2 <- chr_2
out$GL_het_Chr3 <- chr_3
out$GL_het_Chr4 <- chr_4
out$GL_het_Chr5 <- chr_5
out$GL_het_Chr6 <- chr_6
out$GL_het_Chr7 <- chr_7
out$GL_het_Chr8 <- chr_8
out$GL_het_Chr9 <- chr_9
out$GL_het_Chr10 <- chr_10
out$GL_het_Chr11 <- chr_11
out$GL_het_Chr12 <- chr_12
out$GL_het_Chr13 <- chr_13
out$GL_het_Chr14 <- chr_14
out$GL_het_Chr15 <- chr_15
out$GL_het_Chr16 <- chr_16
out$GL_het_Chr17 <- chr_17

write.table(out,"GL_Heterozygosity_dataset.txt",quote=F,col.names=T,row.names=F,sep="\t")
