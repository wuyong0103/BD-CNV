# cnv_bin_gwas_analysis.R

# Usage:
# Rscript cnv_bin_gwas_analysis.R --cnv cnv_file.txt --phe phe_file.txt --chrlen chr_length.txt --binsize 100000 --out results.txt --threads 4

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(parallel))
options(scipen = 999) #Do not use scientific notation

# ----------------------
# Parse command line args
# ----------------------
option_list <- list(
  make_option(c("--cnv"), type="character", help="Input CNV file"),
  make_option(c("--phe"), type="character", help="Phenotype file"),
  make_option(c("--chrlen"), type="character", help="Chromosome length file (2 columns: chr, length)"),
  make_option(c("--binsize"), type="integer", default=100000, help="Bin size [default %default]"),
  make_option(c("--out"), type="character", help="Output file for GWAS results"),
  make_option(c("--threads"), type="integer", default=1, help="Number of parallel threads [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------
# Load input files
# ----------------------
cat("[INFO] Loading files...\n")
cnv <- fread(opt$cnv)
phe <- fread(opt$phe) %>%
  mutate(AFF = as.factor(AFF), SEX = as.factor(SEX), Platform = as.factor(Platform))
chr_sizes <- fread(opt$chrlen, header = FALSE)
colnames(chr_sizes) <- c("CHR", "LENGTH", "V1", "V2")
chr_sizes$CHR <- as.character(chr_sizes$CHR)

# ----------------------
# Generate bins
# ----------------------
cat("[INFO] Generating bins...\n")
bin_size <- opt$binsize
bins <- chr_sizes %>%
  rowwise() %>%
  do({
    chr <- .$CHR
    length <- .$LENGTH
    starts <- seq(1, length, by = bin_size)
    ends <- pmin(starts + bin_size - 1, length)
    data.frame(CHR = chr, BIN_START = starts, BIN_END = ends)
  }) %>%
  ungroup()
bins$BIN_ID <- paste0("chr", bins$CHR, "_", bins$BIN_START, "_", bins$BIN_END)

# ----------------------
# Build sample-bin matrix
# ----------------------
cat("[INFO] Building sample-bin matrix...\n")
samples <- unique(phe$IID)
bin_matrix <- matrix(0, nrow = length(samples), ncol = nrow(bins))
rownames(bin_matrix) <- samples
colnames(bin_matrix) <- bins$BIN_ID

for (i in 1:nrow(cnv)) {
  sid <- cnv$IID[i]
  chr <- as.character(cnv$CHR[i])
  start <- cnv$BP1[i]
  end <- cnv$BP2[i]

  hits <- which(
    bins$CHR == chr &
    !(bins$BIN_END < start | bins$BIN_START > end)
#    bins$CHR == chr & ((start > bins$BIN_START & start < bins$BIN_END) | (end > bins$BIN_START & end < bins$BIN_END))
  )
  
  if (length(hits) > 0 && sid %in% samples) {
    bin_matrix[sid, hits] <- 1
  }
}

# ----------------------
# GWAS by bin (parallel)
# ----------------------
cat("[INFO] Performing bin-level GWAS in parallel...\n")
bin_df <- as.data.frame(bin_matrix)
bin_df$IID <- rownames(bin_df)
merged <- merge(phe, bin_df, by = "IID")
#write.table(merged, "merged.txt", row.names = F)
bin_cols <- grep("^chr", colnames(merged), value=TRUE)

results <- mclapply(bin_cols, function(bin) {
  x <- merged[[bin]]
  if (length(unique(x)) < 2) return(NULL)  # 跳过全 0 或全 1

  # 构建模型公式
  fml <- as.formula(paste0("AFF ~ `", bin, "` + SEX + Platform + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10"))
  model <- tryCatch(glm(fml, data=merged, family=binomial), error=function(e) NULL)

  # GLM 结果
  if (!is.null(model)) {
    coef_summary <- summary(model)$coefficients
    if (bin %in% rownames(coef_summary)) {
      beta <- coef_summary[bin, "Estimate"]
      se <- coef_summary[bin, "Std. Error"]
      p_glm <- coef_summary[bin, "Pr(>|z|)"]
    } else {
      beta <- se <- p_glm <- NA
    }
  } else {
    beta <- se <- p_glm <- NA
  }

  # Fisher exact test
  tab <- table(x, merged$AFF)
  # ensure 2x2
  if (!all(c(0,1) %in% rownames(tab)) || !all(c(1,2) %in% colnames(tab))) {
    new_tab <- matrix(0, nrow=2, ncol=2)
    rownames(new_tab) <- c(0,1)
    colnames(new_tab) <- c(1,2)
    new_tab[rownames(tab), colnames(tab)] <- tab
    tab <- new_tab
  }
  
  a <- tab["1", "2"] + 0.5  # CNV=1 & Case
  b <- tab["1", "1"] + 0.5  # CNV=1 & Control
  c <- tab["0", "2"] + 0.5  # CNV=0 & Case
  d <- tab["0", "1"] + 0.5  # CNV=0 & Control
  a_raw <- tab["1", "2"]
  b_raw <- tab["1", "1"]
  total_case <- sum(merged$AFF == 2)
  total_control <- sum(merged$AFF == 1)

  # OR and SE
  fisher_or <- (a * d) / (b * c)
  fisher_se <- sqrt(1/a + 1/b + 1/c + 1/d)
  fisher_p <- tryCatch(fisher.test(tab)$p.value, error=function(e) NA)

  data.frame(
    BIN = bin,
    Case_CNV = a_raw,
    Total_Case = total_case,
    Ctrl_CNV = b_raw,
    Total_Control = total_control,
    GLM_Beta = beta,
    GLM_SE = se,
    GLM_P = p_glm,
    Fisher_OR = fisher_or,
    Fisher_SE = fisher_se,
    Fisher_P = fisher_p
    )
}, mc.cores=opt$threads)

results <- do.call(rbind, results)


# Save results
# ----------------------
cat("[INFO] Saving results...\n")
cat("[INFO] The number of bins that yield valid regression results:", nrow(results), "\n")
write.table(results, opt$out, sep = "\t", row.names = FALSE, quote = FALSE)

cat("[DONE] Output written to:", opt$out, "\n")
cat("[INFO] Bin GWAS analyzed successfully!\n")
