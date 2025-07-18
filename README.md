# BD CNV study
This Github contains the following two scripts. 
### 1. cnv_bin_gwas_analysis.r
This script is used for conducting genome-wide association study (GWAS) of copy number variant (CNV) based on bins. We divided the chromosomes into consecutive bins, and then used Fisher's test to determine whether the CNVs contained in each bin had a frequency difference between the case and control groups. Odds ratios (ORs), standard errors (SEs), and p-values were calculated for each bin using standard methods for 2×2 tables. Bins with a nominal p-value < 0.05 were considered significant and selected for further investigation. All CNVs falling within significant bins (Fisher’s p-value < 0.05) were aggregated across individuals. CNV regions from different individuals that were separated by ≤1 Mb were merged into unified CNV intervals. <br>
```bash
Rscript cnv_bin_gwas_analysis.r --cnv All_plinked_BDCtrl_All_Del.cnv --phe All_plinked_BDCtrl_All.phe --chrlen hg19_chromosome_length.txt --binsize 1000000 --out Result_Bin_1Mb_del.txt --threads 50
```
The file <I>hg19_chromosome_length.txt</I> contains the length information of each chromosome.
By setting the --binsize parameter, chromosomes can be divided into consecutive bins.
This script supports multi-threading. You can accelerate its execution by setting the --threads parameter.
  
### 2. bin_region_annotation.r
This script is used for annotating CNVs, including annotating CNVs to cytobands as well as genes.

