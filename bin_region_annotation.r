# Usage:
# Rscript Bin_Region_annotation.r --cnv All_plinked_BDCtrl_All_Del --bin Result_Bin_1Mb_del.txt --cyto cytoBand.txt --ann /home/lilab/reference/hg19_gencode/annotation.txt --out Bin_1Mb_del_ranges.txt

# load libraies
suppressMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(optparse)
})

option_list <- list(
  make_option(c("--cnv"), type="character", help="Input CNV file, without postfix"),
  make_option(c("--bin"), type="character", help="Bin GWAS results file"),
  make_option(c("--ann"), type="character", help="Gene annotation file"),
  make_option(c("--out"), type="character", help="Output file for GWAS results"),
  make_option(c("--cyto"), type="character", help="hg19 cytoband file")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- read the data ----
cat("[INFO]Read input...\n")
cnv <- read.table(paste0(opt$cnv, ".cnv"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
fam <- read.table(paste0(opt$cnv, ".fam"), header=FALSE, stringsAsFactors=FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "AFF")
gwas <- read.table(opt$bin, header=TRUE, sep="\t", stringsAsFactors=FALSE)
gene_anno <- read_tsv(opt$ann) %>% filter(!is.na(chr), !is.na(start), !is.na(end), !is.na(gene_type), !is.na(gene_name), !is.na(ensembl))
cyto <- read.table(opt$cyto, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(cyto) <- c("chr", "start", "end", "band", "stain")

# ---- merge fam info to CNV ----
cnv <- merge(cnv, fam[, c("FID", "IID", "AFF")], by=c("FID", "IID"))

# ---- set the total number of sample size ----
total_case <- 4842
total_control <- 7820

# ---- extract the significant bin and its position of chromosome ----
sig_bins <- subset(gwas, Fisher_P < 0.05)
bin_info <- do.call(rbind, strsplit(sig_bins$BIN, "_"))
colnames(bin_info) <- c("CHR_STR", "START", "END")
sig_bins$CHR <- as.integer(gsub("chr", "", bin_info[,1]))
sig_bins$BIN_START <- as.integer(bin_info[,2])
sig_bins$BIN_END <- as.integer(bin_info[,3])

# ---- Extract CNVs faling in the significant bins ----
all_cnv_sig_bin <- data.frame()
for (i in seq_len(nrow(sig_bins))) {
  bin <- sig_bins[i, ]
  chr <- bin$CHR
  start <- bin$BIN_START
  end <- bin$BIN_END

  cnv_in_bin <- cnv %>%
    filter(CHR == chr & BP1 <= end & BP2 >= start)
  all_cnv_sig_bin <- rbind(all_cnv_sig_bin, cnv_in_bin)
}
all_cnv_sig_bin <- unique(all_cnv_sig_bin)

# ---- Merge the overlap CNV regions ----
cat("[INFO]Build cnv GRanges...\n")
all_cnv_sig_bin <- all_cnv_sig_bin %>%
  filter(!is.na(CHR), !is.na(BP1), !is.na(BP2), !is.na(IID), !is.na(AFF))
cnv_gr <- GRanges(seqnames = paste0("chr", all_cnv_sig_bin$CHR),
                  ranges = IRanges(start = all_cnv_sig_bin$BP1, end = all_cnv_sig_bin$BP2),
                  sample_id = all_cnv_sig_bin$IID,
                  aff = all_cnv_sig_bin$AFF)

cat("[INFO]Merge adjacent CNV with distance < 1Mb...\n")
merged_gr <- reduce(cnv_gr, min.gapwidth = 1e6)

cat("[INFO]Build cytoband GRanges...\n")
# ---- build cytoband GRanges ----
cyto_gr <- GRanges(
  seqnames = cyto$chr,
  ranges = IRanges(start = cyto$start, end = cyto$end),
  band = cyto$band
)

# ---- Building mapping ----
cat("[INFO]Build mapping...\n")
results <- list()

for (i in seq_along(merged_gr)) {
  region <- merged_gr[i]
  chr <- gsub("chr", "", as.character(seqnames(region)))
  region_start <- start(region)
  region_end <- end(region)

  # Find overlap CNV
  overlaps <- findOverlaps(region, cnv_gr)
  cnv_in_region <- cnv_gr[subjectHits(overlaps)]

  # Extract samples
  df <- as.data.frame(cnv_in_region)
  case_ids <- unique(df$sample_id[df$aff == 2])
  ctrl_ids <- unique(df$sample_id[df$aff == 1])

  a <- length(case_ids)
  b <- length(ctrl_ids)
  c <- total_case - a
  d <- total_control - b

  # Fisher test
  or <- NA
  if (all(!is.na(c(a,b,c,d))) && all(c(a,b,c,d) >= 0)) {
    or <- if (a > 0 && b > 0 && c > 0 && d > 0) (a * d)/(b * c) else NA
    pval <- fisher.test(matrix(c(a, b, c, d), nrow=2))$p.value
  } else {
    pval <- NA
  }

  # ---- annotate genes ----
  gene_gr <- GRanges(seqnames = paste0("chr", gene_anno$chr),
                     ranges = IRanges(start = gene_anno$start, end = gene_anno$end),
                     gene_name = gene_anno$gene_name,
                     gene_type = gene_anno$gene_type,
                     ensembl_id = gene_anno$ensembl)

  gene_ol <- findOverlaps(region, gene_gr)
  if (length(gene_ol) > 0) {
    gdf <- as.data.frame(gene_gr[subjectHits(gene_ol)])
    genes <- paste(unique(gdf$gene_name), collapse = ";")
    gtypes <- paste(unique(gdf$gene_type), collapse = ";")
    gids <- paste(unique(gdf$ensembl_id), collapse = ";")
  } else {
    genes <- NA
    gtypes <- NA
    gids <- NA
  }

  # ------ cytoband annotation --------
  cyto_ol <- findOverlaps(region, cyto_gr)
  if (length(cyto_ol) > 0) {
    band_names <- unique(mcols(cyto_gr[subjectHits(cyto_ol)])$band)
    bands <- paste(band_names, collapse = ";")
  } else {
    bands <- NA
  }

  results[[length(results)+1]] <- data.frame(
    CHR = chr,
    CNV_START = region_start,
    CNV_END = region_end,
    Cytoband = paste0(chr, bands),
    Case_with_CNV = a,
    Control_with_CNV = b,
    Fisher_OR = or,
    Fisher_P = pval,
    Genes = genes,
#    Gene_Types = gtypes,
    Ensembl_IDs = gids
  )
}

# ---- output ----
cat("[INFO]Write to output...\n")
final_result <- do.call(rbind, results)
write.table(final_result, opt$out, sep="\t", quote=FALSE, row.names=FALSE)
cat("[Done]Finsihed successfully!\n")
