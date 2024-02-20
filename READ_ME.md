This repository contains the code relative to the RADseq analysis (and part of the microsatellites analysis) for the Northern Elephant Seal paper from Hoffman et al. (2024)
The repo contains also the code to run SLiM simulations on genetic load.

1. Bioinformatics.sh describe the procedure to get form raw fastq sequencing data to quality filtered SNP genotypes;
2. sMLH_g2.r describes how to get sMLH and g2 values for SNPs and microsatellites (include also code for calculating
   sMLH separately for each chromosome and within sliding windows);
3. GenLikelihood_Heterozygosity.sh descirbes how to get genotype likelihood based heterozigosity for each sample.
4. PCA_script.r is to run PCA on SNP and microsatellites data;
5. The 'Demography' folder contains fastsimcoal2 .tpl and .est files used for the three tested demographic scenarios
   and the workflow ('Demography.sh') to run fastsimcoal2, select the best runs and generate 95% CI based on non-parametric bootstrap;
6. The 'SLiM' folder contains everything necessary to re-run the SLiM simulations described in the paper. More details are present in the READ_ME files
   in the folder (WF = Wright-Fisher models).
   