# BD CNV study
This Github contains only one script. 
### 1. cnv_gene_gwas.R
This script is used for conducting genome-wide association study (GWAS) of copy number variant (CNV) based on genes. A gene was considered affected in an individual if at least one CNV overlapped its genomic region with at least one base pair. We then constructed a binary matrix, with rows representing individuals and columns representing genes, where each entry indicated the presence (1) or absence (0) of CNVs affecting the corresponding gene in each individual. For each gene, we performed logistic regression to test for association between CNV burden and case-control status, adjusting for sex, genotyping platform, and the top 10 principal components derived from SNP genotype data to account for population stratification. In cases of sparse data or complete separation, Firth’s penalized logistic regression was applied using the <I>logistf</I> package in R to obtain bias-reduced estimates. Genes with at least one CNV-affected individuals were included in the analysis. <br>
```bash
Rscript cnv_gene_gwas.R \
    --phe Phe.txt \
    --cnv CNV_Plink_Format_Del.cnv \
    --gene GeneAnnotation.txt \
    --out CNV_GWAS_Result_Del.txt \
    --threads 1
```
The `Phe.txt` has the following columns: <br>
```
FID     IID     AFF     SEX     Platform        C1      C2      C3      C4      C5      C6      C7      C8      C9      C10
Sample1      Sample1      1       2       ASA     -0.0022 -0.0088 -0.0004 -0.0048 0.0029  0.001   0.0023  -0.0028 0.0024  -0.0117
Sample2      Sample2      1       1       ASA     0       -0.0027 -0.0001 -0.0041 -0.0064 0.0071  -0.0014 0.0069  0.0029  -0.0086
Sample3      Sample3      1       1       ASA     -0.0015 -0.0064 -0.0021 -0.0117 0.0006  -0.0017 0.002   -0.0062 -0.0056 -0.0139
```
Case: AFF=1; Control: AFF=0

The `Phe.txt` has the following columns: <br>
```
10      97737121        97760907        ENSG00000155256
10      97766751        97771999        ENSG00000120057
10      97849843        97871580        ENSG00000155265
```

This script supports multi-threading. You can accelerate its execution by setting the `--threads` parameter (Note, this option consumes an excessive amount of memory. It is suggested to split the GeneAnnotation.txt file, run cnv_gene_gwas.R separately, and then merge the results.). <br>
