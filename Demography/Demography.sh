## Remove from .bam files sex chromosome scaffolds 

for i in *.bam
do
bedtools intersect -a $i -b grep_Sex.bed -ubam -v > ${i%.bam}_NoSex.bam
done

## Create SFS with angsd + 100 bootstrap
#!/bin/bash
# Do saf
angsd -uniqueOnly 1 -remove_bads 1 -baq 1 -minMapQ 20 -minQ 30 -minInd 96 -setMinDepth 288 -setMaxDepth 1600 -doCounts 1 -GL 1 -out NES_angsd -anc Mirounga_angustirostris_HiC.fasta -ref Mirounga_angustirostris_HiC.fasta -bam bam_files.txt -doSaf 1 -P 8

# Do folded SFS
realSFS NES_angsd.saf.idx -maxIter 100 -fold 1 -P 8 > NES_folded_SFS.sfs

# Bootstrap 100
realSFS NES_angsd.saf.idx -maxIter 100 -bootstrap 100 -fold 1 -P 8 > NES_folded_BootstrappedSFS.sfs

## Run the models based on the empirical full SFS: 
# Null
for i in {1..100};
do
fsc2702 -t Null.tpl -n 100000 -m -e Null.est -M -L 40 -q -w 0.01 --foldedSFS -x -C 5 --nosingleton -c 4
mv Null Null_$i
done

# Bott_6
for i in {1..100};
do
fsc2702 -t Bott_6.tpl -n 100000 -m -e Bott_6.est -M -L 40 -q -w 0.01 --foldedSFS -x -C 5 --nosingleton -c 4
mv Bott_6 Bott_6_$i
done

# Bott_10
for i in {1..100};
do
fsc2702 -t Bott_10.tpl -n 100000 -m -e Bott_10.est -M -L 40 -q -w 0.01 --foldedSFS -x -C 5 --nosingleton -c 4
mv Bott_10 Bott_10_$i
done

# Non-parametric bootstrap (from ANGSD bootstrapped SFS):
# First format ANGSD bootstrap output and place each bootstrap replicate in a separate folder:
for i in {1..100}
do
head -$i NES_folded_BootstrappedSFS_900more.sfs | tail -1 > Boot.sfs
cat header.txt Boot.sfs > Bootstrap/Boot_$i/Null_MAFpop0.obs # Change name depending on the model so that it math=ches with .tpl and .est files!
rm Boot.sfs
done
# Where header.txt is:
# 1 observations
# d0_0 d0_1 d0_2 d0_3 d0_4 d0_5 d0_6 d0_7 d0_8 d0_9 d0_10 d0_11 d0_12 d0_13 d0_14 d0_15 d0_16 d0_17 d0_18 d0_19 d0_20 d0_21 d0_22 d0_23 d0_24 d0_25 d0_26 d0_27 d0_28 d0_29 d0_30 d0_31 d0_32 d0_33 d0_34 d0_35 d0_36 d0_37 d0_38 d0_39 d0_40 d0_41 d0_42 d0_43 d0_44 d0_45 d0_46 d0_47 d0_48 d0_49 d0_50 d0_51 d0_52 d0_53 d0_54 d0_55 d0_56 d0_57 d0_58 d0_59 d0_60 d0_61 d0_62 d0_63 d0_64 d0_65 d0_66 d0_67 d0_68 d0_69 d0_70 d0_71 d0_72 d0_73 d0_74 d0_75 d0_76 d0_77 d0_78 d0_79 d0_80 d0_81 d0_82 d0_83 d0_84 d0_85 d0_86 d0_87 d0_88 d0_89 d0_90 d0_91 d0_92 d0_93 d0_94 d0_95 d0_96 d0_97 d0_98 d0_99 d0_100 d0_101 d0_102 d0_103 d0_104 d0_105 d0_106 d0_107 d0_108 d0_109 d0_110 d0_111 d0_112 d0_113 d0_114 d0_115 d0_116 d0_117 d0_118 d0_119 d0_120 d0_121 d0_122 d0_123 d0_124 d0_125 d0_126 d0_127 d0_128 d0_129 d0_130 d0_131 d0_132 d0_133 d0_134 d0_135 d0_136 d0_137 d0_138 d0_139 d0_140 d0_141 d0_142 d0_143 d0_144 d0_145 d0_146 d0_147 d0_148 d0_149 d0_150 d0_151 d0_152 d0_153 d0_154 d0_155 d0_156 d0_157 d0_158 d0_159 d0_160 d0_161 d0_162 d0_163 d0_164 d0_165 d0_166 d0_167 d0_168 d0_169 d0_170 d0_171 d0_172 d0_173 d0_174 d0_175 d0_176 d0_177 d0_178 d0_179 d0_180 d0_181 d0_182 d0_183 d0_184 d0_185 d0_186 d0_187 d0_188 d0_189 d0_190 d0_191 d0_192

# Then re-run fsc.
# Null 
for j in {1..100};
do
cd Boot_$j

for i in {1..100};
do
fsc2702 -t Null.tpl -n 100000 -m -e Null.est -M -L 40 -q -w 0.01 --foldedSFS -x -C 5 --nosingleton -c 6
mv Null Null_$i
done

cd ../
done

# Bott_6 
for j in {1..100};
do
cd Boot_$j

for i in {1..100};
do
fsc2702 -t Bott_6.tpl -n 100000 -m -e Bott_6.est -M -L 40 -q -w 0.01 --foldedSFS -x -C 5 --nosingleton -c 6
mv Bott_6 Bott_6_$i
done

cd ../
done

# Bott_10 
for j in {1..100};
do
cd Boot_$j

for i in {1..100};
do
fsc2702 -t Bott_10.tpl -n 100000 -m -e Bott_10.est -M -L 40 -q -w 0.01 --foldedSFS -x -C 5 --nosingleton -c 6
mv Bott_10 Bott_10_$i
done

cd ../
done

# Get best runs from each bootstrap: 
for j in {1..100}
do

cd Boot_$j

head -1 LGM_1_1/LGM_1.bestlhoods > Header.txt

for i in {1..20}
do
tail -1 LGM_1_$i/LGM_1.bestlhoods >> Lhoods.txt
done

cat Header.txt Lhoods.txt > BestLhoods.txt
rm Lhoods.txt Header.txt
 
Rscript get_best.r
## get_best.r
#data <- read.table("BestLhoods.txt", h = T)
#ind <- which(data$MaxEstLhood == max(data$MaxEstLhood))
#out <- data[ind,]
#write.table(out, "bestL.txt", quote = F, col.names = F, row.names = F, sep = '\t')

cd ../
done

# Merge them into a single file:
cat Header.txt Boot_*/bestL.txt > Non_Param_boots.txt

# Get 95% CIs:
data<-read.table("Non_Param_boots.txt",h=T)

for (i in c(1:3,5)){
data[,i]<-data[,i]/2
}

apply(data,2,function(x) quantile(x, probs=c(0.025,0.975)))

