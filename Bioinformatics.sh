### From raw reads to SNPs

# Run fastqc to check reads quality: need trimming!

# map to genome
for i in *.1.fq; do ~/bin/bwa-0.7.17/bwa mem -t 5 M_angustirostris.fasta $i ${i%.1.fq}.2.fq > ${i%.1.fq}.sam; done

# process sam files
for i in *.sam; do samtools view -SbF 4 $i > ${i%sam}filt.bam; done
for i in *.filt.bam; do gatk SortSam -I $i -O ${i%bam}sorted.bam -SO coordinate; done
for i in *.sorted.bam; do gatk AddOrReplaceReadGroups -I $i -O ${i%bam}RG.bam -ID ${i%.bam} -PL Illumina -LB RADseq -PU HWI-1KL168:100:HNJT5BCXX -SM ${i%.filt.sorted.bam} --VALIDATION_STRINGENCY SILENT; done
for i in *.RG.bam; do gatk MarkDuplicates --REMOVE_DUPLICATES true --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT -I $i -O ${i%bam}rmdup.bam -M ${i%bam}rmdupmetrics; done

# Merge everything into a single bam file (with gatk MergeSamFiles) and index it (samtools index)

# Call SNPs with gatk
# first create a bunch of required files
samtools faidx M_angustirostris.fasta
gatk CreateSequenceDictionary -R M_angustirostris.fasta -O M_angustirostris.dict

gatk HaplotypeCaller -I Merged_bams.bam -R M_angustirostris.fasta -O NES_rawSNPs.vcf -stand_call_conf 30

# Filter SNPs (418,832)
vcftools --vcf NES_rawSNPs.vcf --min-alleles 2 --max-alleles 2 --remove-indels --minDP 5 --minGQ 5 --recode --out NES_DP5GQ5 #321,124
vcftools --vcf NES_DP5GQ5.recode.vcf --max-missing 0.5 --out NES_DP5GQ5_MD50 #41,605
vcftools --vcf NES_DP5GQ5_MD50.recode.vcf --max-meanDP 19.06 --recode --out NES_DP5GQ5_MD50_maxDP #38,668
vcftools --vcf NES_DP5GQ5_MD50_maxDP.recode.vcf --maf 0.01 --recode --out NES_DP5GQ5_MD50_maxDP_MAF01 #15,229
# Remove samples with more than 50% missing genotypes.
vcftools --vcf NES_DP5GQ5_MD50_maxDP_MAF01.recode.vcf --missing-indv --out NES_Missingness
# Switch to R
R
data<-read.table('NES_Missingness.imiss',h=T)
ind<-which(data$Miss<(15229*0.5))
out<-data[ind,1]
write.table(out,'SamplesToRemove.txt',quote=F,col.names=F,row.names=F,sep='\t')
q()
n
# Back to bash
vcftools --vcf NES_DP5GQ5_MD50_maxDP_MAF01.recode.vcf --remove SamplesToRemove.txt --recode --out NES_Final
# Convert to plink
vcftools --NES_Final.recode.vcf --plink --out NES_Final
plink --file NES_Final --recodeA --aec --out NES_Final


