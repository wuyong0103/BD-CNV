#!/usr/bin/env Rscript

# ==== 载入库 ====
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(epitools)
  library(optparse)
  library(parallel)
  library(logistf)
})

# ==== 命令行参数 ====
option_list <- list(
  make_option(c("-p", "--phe"), type = "character", help = "Phenotype file (e.g. phe.txt)", metavar = "FILE"),
  make_option(c("-c", "--cnv"), type = "character", help = "CNV segments file (e.g. cnv_segments.txt)", metavar = "FILE"),
  make_option(c("-g", "--gene"), type = "character", help = "Gene annotation file", metavar = "FILE"),
  make_option(c("-o", "--out"), type = "character", help = "Output result file", metavar = "FILE"),
  make_option(c("-t", "--threads"), type = "integer", default = floor(parallel::detectCores() / 2),
              help = "Number of threads to use [default: %default]", metavar = "INT")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("Using", opt$threads, "threads\n")

# ==== Step 1: 读取数据 ====
phe <- fread(opt$phe)
cnv <- fread(opt$cnv)
genes <- fread(opt$gene, col.names = c("CHR", "START", "END", "GENE"))

# ==== Step 2: 扩展基因窗口 ====
genes <- genes %>%
  mutate(WIN_START = pmax(START, 0),
         WIN_END = END)

# ==== Step 3: CNV区段格式处理 ====
cnv <- cnv %>%
  rename(CNV_START = BP1, CNV_END = BP2)

setDT(cnv)
setDT(genes)
setkey(genes, CHR, WIN_START, WIN_END)

# ==== Step 4: 区间重叠 ====
overlap <- foverlaps(cnv, genes, by.x = c("CHR", "CNV_START", "CNV_END"), type = "any", nomatch = NULL)
gene_cnv_presence <- overlap[, .(CNV_PRESENT = 1), by = .(IID, GENE)]

# ==== Step 5: 所有样本 × 所有基因 ====
all_pairs <- expand.grid(IID = phe$IID, GENE = unique(genes$GENE), stringsAsFactors = FALSE)
gene_matrix <- merge(all_pairs, gene_cnv_presence, by = c("IID", "GENE"), all.x = TRUE)
gene_matrix$CNV_PRESENT[is.na(gene_matrix$CNV_PRESENT)] <- 0

# ==== Step 6: 合并表型 ====
gene_matrix <- merge(gene_matrix, phe, by = "IID")

# ==== Step 7: 基因分析函数 ====
analyze_gene <- function(gene_id, df) {
  cat(paste0("[INFO] processing ", gene_id, " ...\n"))
  sub <- df[df$GENE == gene_id, ]
  print(table(sub$CNV_PRESENT))
  if (length(unique(sub$CNV_PRESENT)) < 2) return(NULL)
  
  tab <- table(sub$AFF, sub$CNV_PRESENT)
  print(tab)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  
  fisher <- fisher.test(tab)
  fisher_OR <- as.numeric(fisher$estimate)
  fisher_CI <- fisher$conf.int
  fisher_p <- fisher$p.value

  sub$SEX <- as.factor(sub$SEX)
  sub$Platform <- as.factor(sub$Platform)

  if (any(tab == 0)) {
    glm_fit <- try(logistf(AFF ~ CNV_PRESENT + SEX + Platform + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10, data = sub), silent = TRUE)
    firth_model <- TRUE
  } else {
    glm_fit <- try(glm(AFF ~ CNV_PRESENT + SEX + Platform + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10, data = sub, family = binomial()), silent = TRUE)
    firth_model <- FALSE
  }

  if (inherits(glm_fit, "try-error")) return(NULL)

  if (firth_model) {
    OR_logit <- exp(coef(glm_fit)["CNV_PRESENT"])
    confint_logit <- exp(confint(glm_fit)["CNV_PRESENT", ])
    p_logit <- glm_fit$prob["CNV_PRESENT"]
  } else {
    OR_logit <- exp(coef(glm_fit)["CNV_PRESENT"])
    confint_logit <- tryCatch({
      exp(confint(glm_fit)["CNV_PRESENT", ])
    }, error = function(e) c(NA, NA))
    p_logit <- summary(glm_fit)$coefficients["CNV_PRESENT", "Pr(>|z|)"]
  }

  num_case <- sum(sub$CNV_PRESENT == 1 & sub$AFF == 1)
  num_control <- sum(sub$CNV_PRESENT == 1 & sub$AFF == 0)
      
  data.frame(
    Gene = gene_id,
    N_case_with_CNV = num_case,
    N_control_with_CNV = num_control,
    Fisher_OR = fisher_OR,
    Fisher_LCI = fisher_CI[1],
    Fisher_UCI = fisher_CI[2],
    Fisher_p = fisher_p,
    Logit_OR = OR_logit,
    Logit_LCI = confint_logit[1],
    Logit_UCI = confint_logit[2],
    Logit_p = p_logit
  )
}

# ==== Step 8: 并行分析所有基因 ====
gene_list <- unique(gene_matrix$GENE)
results_list <- mclapply(gene_list, function(g) analyze_gene(g, gene_matrix), mc.cores = opt$threads)
#results_list <- lapply(gene_list, function(g) analyze_gene(g, gene_matrix))
results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

# ==== Step 9: 输出 ====
results <- results[order(Fisher_p)]
fwrite(results, file = opt$out, sep = "\t", quote = FALSE)
