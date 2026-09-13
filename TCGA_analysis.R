###############################################################################
# TCGA_analysis.R
#
# TCGA-LUSC transcriptomic biological relevance and survival analyses for:
# "Prognostic Implications and Proteomic Classification of
#  Keratinizing versus Non-Keratinizing Lung Squamous Cell Carcinoma"
#
# Required input files:
#   1. gencode.v36.annotation.gtf.gene.probemap
#   2. TCGA-LUSC.star_tpm.tsv.gz
#   3. TCGA-LUSC.survival.tsv
#
# This script uses the finalized clean 10-protein logistic-model coefficients
# from the proteomic discovery cohort. Because TCGA RNA-seq and LFQ protein
# abundance are measured on different scales, the TCGA-derived KSCC score is
# treated as an ordinal molecular index rather than a calibrated probability.
#
# Main analyses:
#   - TCGA-derived KSCC score
#   - Spearman correlation with TP63
#   - Spearman correlation with 8 keratinization-related KRT genes
#   - Benjamini-Hochberg correction across all 9 correlation tests
#   - Kaplan-Meier analysis after median-based score stratification
#   - Cox model for KSCC score per 1-SD increase
###############################################################################


###############################################################################
# 0. Packages
###############################################################################

library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(survival)
library(survminer)


###############################################################################
# 1. Read TCGA-LUSC TPM data and gene annotation
###############################################################################

probemap <- read.table(
  "gencode.v36.annotation.gtf.gene.probemap",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

tpm_data <- read.table(
  "TCGA-LUSC.star_tpm.tsv.gz",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  all(c("id", "gene") %in% colnames(probemap))
)

final_tpm <- merge(
  probemap[, c("id", "gene")],
  tpm_data,
  by.x = "id",
  by.y = colnames(tpm_data)[1]
)


###############################################################################
# 2. Final clean 10-protein logistic-model coefficients
#
# These are the coefficients from the finalized 97-sample clean model:
# Group ~ ERGIC2 + PKP1 + PRXL2B + IDO1 + LGALS7B +
#         SIGLEC1 + COL5A1 + MRPL2 + GAMT + CPNE2
#
# KSCC is the modeled event.
###############################################################################

intercept <- 27.40185333

coeffs <- c(
  ERGIC2  =  0.07064915,
  PKP1    = -0.15902274,
  PRXL2B  = -0.46972149,
  IDO1    = -0.67017050,
  LGALS7B =  0.13301533,
  SIGLEC1 = -0.16987744,
  COL5A1  =  0.18125215,
  MRPL2   = -0.51363358,
  GAMT    =  0.67294613,
  CPNE2   = -0.69540355
)

model_genes <- names(coeffs)

stopifnot(
  length(model_genes) == 10,
  identical(
    model_genes,
    c(
      "ERGIC2", "PKP1", "PRXL2B", "IDO1", "LGALS7B",
      "SIGLEC1", "COL5A1", "MRPL2", "GAMT", "CPNE2"
    )
  )
)


###############################################################################
# 3. Gene mapping
#
# If LGALS7B is unavailable in the TCGA transcriptomic annotation but LGALS7
# is present, LGALS7 is used as the transcriptomic proxy, matching the
# finalized analysis.
###############################################################################

gene_map <- setNames(
  model_genes,
  model_genes
)

if (
  "LGALS7B" %in% model_genes &&
  !("LGALS7B" %in% final_tpm$gene) &&
  ("LGALS7" %in% final_tpm$gene)
) {
  gene_map["LGALS7B"] <- "LGALS7"
}

krt_genes <- c(
  "KRT5",
  "KRT6A",
  "KRT6B",
  "KRT13",
  "KRT14",
  "KRT15",
  "KRT16",
  "KRT17"
)

mechanism_genes <- c(
  "TP63",
  krt_genes
)

target_genes_tcga <- unique(
  c(
    unname(gene_map),
    mechanism_genes
  )
)

missing_model_genes <- setdiff(
  unname(gene_map),
  final_tpm$gene
)

if (length(missing_model_genes) > 0) {
  stop(
    paste(
      "These model genes are missing in TCGA:",
      paste(missing_model_genes, collapse = ", ")
    )
  )
}

missing_mechanism_genes <- setdiff(
  mechanism_genes,
  final_tpm$gene
)

if (length(missing_mechanism_genes) > 0) {
  stop(
    paste(
      "These TP63/KRT genes are missing in TCGA:",
      paste(missing_mechanism_genes, collapse = ", ")
    )
  )
}

cat("\nModel gene mapping used for TCGA:\n")
print(gene_map)


###############################################################################
# 4. Build TCGA expression matrix and retain primary tumors only
###############################################################################

tcga_expr <- final_tpm %>%
  filter(gene %in% target_genes_tcga) %>%
  group_by(gene) %>%
  summarise(
    across(where(is.numeric), mean),
    .groups = "drop"
  ) %>%
  column_to_rownames("gene")

tcga_expr_t <- as.data.frame(
  t(tcga_expr),
  check.names = FALSE
)

# TCGA sample-type code 01 = primary tumor
tcga_expr_t$sample_type <- substr(
  rownames(tcga_expr_t),
  14,
  15
)

tcga_tumor <- tcga_expr_t %>%
  filter(sample_type == "01")

tcga_tumor$sample_type <- NULL

cat(
  "\nTCGA primary tumor samples retained:",
  nrow(tcga_tumor),
  "\n"
)

# The finalized analysis used 501 TCGA-LUSC primary tumor samples
stopifnot(
  nrow(tcga_tumor) == 501
)

tcga_tumor_log2 <- log2(
  tcga_tumor + 1
)


###############################################################################
# 5. Construct TCGA model-gene matrix and calculate KSCC score
#
# The 10 transcriptomic features are standardized within TCGA-LUSC before the
# proteomic-model coefficients are applied. The resulting linear predictor is
# treated as an ordinal molecular index, not a calibrated probability.
###############################################################################

tcga_model_expr <- data.frame(
  row.names = rownames(tcga_tumor_log2)
)

for (model_gene in model_genes) {
  tcga_gene <- unname(
    gene_map[model_gene]
  )

  tcga_model_expr[[model_gene]] <-
    tcga_tumor_log2[[tcga_gene]]
}

if (anyNA(tcga_model_expr)) {
  stop("Missing values detected in TCGA model-gene expression matrix.")
}

tcga_model_expr_scaled <- as.data.frame(
  scale(tcga_model_expr)
)

if (anyNA(tcga_model_expr_scaled)) {
  stop("NA generated during TCGA Z-score scaling. Check gene variance.")
}

tcga_model_expr_scaled$KSCC_Score <-
  intercept +
  as.numeric(
    as.matrix(
      tcga_model_expr_scaled[
        ,
        model_genes,
        drop = FALSE
      ]
    ) %*%
      as.numeric(
        coeffs[model_genes]
      )
  )

tcga_score_df <- tcga_model_expr_scaled

cat("\nTCGA KSCC score summary:\n")
print(
  summary(
    tcga_score_df$KSCC_Score
  )
)


###############################################################################
# 6. Spearman correlations:
#    KSCC score vs TP63 and 8 keratinization-related KRT genes
#
# BH correction is applied jointly across all 9 prespecified tests.
###############################################################################

correlation_genes <- c(
  "TP63",
  krt_genes
)

correlation_results <- lapply(
  correlation_genes,
  function(gene) {

    res <- cor.test(
      tcga_score_df$KSCC_Score,
      tcga_tumor_log2[
        rownames(tcga_score_df),
        gene
      ],
      method = "spearman",
      exact = FALSE
    )

    data.frame(
      Gene = gene,
      Spearman_rho = as.numeric(res$estimate),
      P_value = res$p.value,
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

correlation_results$Adjusted_P_value <- p.adjust(
  correlation_results$P_value,
  method = "BH"
)

correlation_results$Significance <- case_when(
  correlation_results$Adjusted_P_value < 0.001 ~ "***",
  correlation_results$Adjusted_P_value < 0.01 ~ "**",
  correlation_results$Adjusted_P_value < 0.05 ~ "*",
  TRUE ~ "NS"
)

cat(
  "\n============================================================\n"
)
cat("TCGA CORRELATION RESULTS\n")
cat(
  "============================================================\n"
)

print(
  correlation_results
)

write.csv(
  correlation_results,
  "TCGA_KSCC_Score_Correlation_Results.csv",
  row.names = FALSE
)


###############################################################################
# 7. TP63 scatter plot
###############################################################################

plot_tp63_df <- data.frame(
  KSCC_Score =
    tcga_score_df$KSCC_Score,
  TP63 =
    tcga_tumor_log2[
      rownames(tcga_score_df),
      "TP63"
    ]
)

p_tp63 <- ggplot(
  plot_tp63_df,
  aes(
    x = KSCC_Score,
    y = TP63
  )
) +
  geom_point(
    color = "grey40",
    alpha = 0.5,
    size = 1.8
  ) +
  geom_smooth(
    method = "lm",
    color = "#D62728",
    fill = "#D62728",
    alpha = 0.15,
    linewidth = 1.2
  ) +
  stat_cor(
    method = "spearman",
    label.x.npc = "right",
    label.y.npc = "bottom",
    size = 5,
    fontface = "bold"
  ) +
  theme_classic(
    base_size = 14
  ) +
  labs(
    title =
      "Association between KSCC score and TP63",
    x =
      "TCGA-derived KSCC score",
    y =
      "TP63 expression (log2 TPM)"
  )

ggsave(
  "Figure_TCGA_KSCC_Score_TP63.pdf",
  p_tp63,
  width = 5,
  height = 5
)


###############################################################################
# 8. KRT correlation heatmap
###############################################################################

krt_results <- correlation_results %>%
  filter(Gene %in% krt_genes) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = krt_genes
    )
  ) %>%
  arrange(Gene)

cor_matrix_vertical <- matrix(
  krt_results$Spearman_rho,
  ncol = 1
)

rownames(
  cor_matrix_vertical
) <- as.character(
  krt_results$Gene
)

colnames(
  cor_matrix_vertical
) <- "KSCC score"

label_matrix_vertical <- matrix(
  paste0(
    round(
      krt_results$Spearman_rho,
      3
    ),
    "\n",
    krt_results$Significance
  ),
  ncol = 1
)

rownames(
  label_matrix_vertical
) <- as.character(
  krt_results$Gene
)

colnames(
  label_matrix_vertical
) <- "KSCC score"

positive_red_colors <- colorRampPalette(
  c(
    "white",
    "#F4A6A6",
    "#D62728",
    "#7F0000"
  )
)(100)

pheatmap(
  cor_matrix_vertical,
  color = positive_red_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = label_matrix_vertical,
  number_color = "black",
  fontsize_number = 10,
  fontsize_row = 11,
  fontsize_col = 12,
  angle_col = 0,
  border_color = "white",
  main =
    "KSCC score and keratinization-related genes",
  filename =
    "Figure_TCGA_KSCC_Score_KRT_Association_Heatmap.pdf",
  width = 3.8,
  height = 5.5
)


###############################################################################
# 9. Save TCGA sample-level KSCC scores
###############################################################################

tcga_output <- data.frame(
  Sample =
    rownames(tcga_score_df),
  Patient =
    substr(
      rownames(tcga_score_df),
      1,
      12
    ),
  KSCC_Score =
    tcga_score_df$KSCC_Score,
  stringsAsFactors = FALSE
)

write.csv(
  tcga_output,
  "TCGA_LUSC_KSCC_Score.csv",
  row.names = FALSE
)


###############################################################################
# 10. TCGA-LUSC overall survival data
###############################################################################

score_patient <- tcga_output %>%
  group_by(Patient) %>%
  summarise(
    KSCC_Score =
      mean(
        KSCC_Score,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

surv_raw <- read.table(
  "TCGA-LUSC.survival.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_survival_cols <- c(
  "sample",
  "OS",
  "OS.time"
)

missing_survival_cols <- setdiff(
  required_survival_cols,
  colnames(surv_raw)
)

if (length(missing_survival_cols) > 0) {
  stop(
    paste(
      "These expected survival columns are missing:",
      paste(
        missing_survival_cols,
        collapse = ", "
      )
    )
  )
}

surv_clean <- surv_raw %>%
  transmute(
    Patient =
      substr(
        as.character(sample),
        1,
        12
      ),
    OS_status =
      suppressWarnings(
        as.numeric(OS)
      ),
    OS_time =
      suppressWarnings(
        as.numeric(OS.time)
      )
  ) %>%
  filter(
    !is.na(Patient),
    !is.na(OS_status),
    !is.na(OS_time),
    OS_time > 0,
    OS_status %in% c(0, 1)
  ) %>%
  distinct(
    Patient,
    .keep_all = TRUE
  )

surv_df <- inner_join(
  score_patient,
  surv_clean,
  by = "Patient"
)

cat(
  "\n============================================================\n"
)
cat("FINAL MATCHED TCGA-LUSC SURVIVAL DATASET\n")
cat(
  "============================================================\n"
)

cat(
  "N =",
  nrow(surv_df),
  "\n"
)

cat(
  "Deaths =",
  sum(
    surv_df$OS_status == 1
  ),
  "\n"
)

cat(
  "Censored =",
  sum(
    surv_df$OS_status == 0
  ),
  "\n"
)

# Finalized manuscript/reviewer-response cohort
stopifnot(
  nrow(surv_df) == 493,
  sum(surv_df$OS_status == 1) == 211
)


###############################################################################
# 11. Continuous KSCC score and overall survival
#     HR is reported per 1-SD increase.
###############################################################################

surv_df$KSCC_Score_z <- as.numeric(
  scale(
    surv_df$KSCC_Score
  )
)

cox_cont <- coxph(
  Surv(
    OS_time,
    OS_status
  ) ~ KSCC_Score_z,
  data = surv_df
)

cox_cont_sum <- summary(
  cox_cont
)

HR_cont <- cox_cont_sum$coefficients[
  1,
  "exp(coef)"
]

P_cont <- cox_cont_sum$coefficients[
  1,
  "Pr(>|z|)"
]

CI_cont <- exp(
  confint(
    cox_cont
  )
)

ph_cont <- cox.zph(
  cox_cont
)

cat(
  "\n============================================================\n"
)
cat("CONTINUOUS KSCC SCORE\n")
cat(
  "============================================================\n"
)

cat(
  sprintf(
    "HR per 1-SD increase = %.3f\n",
    HR_cont
  )
)

cat(
  sprintf(
    "95%% CI = %.3f - %.3f\n",
    CI_cont[1],
    CI_cont[2]
  )
)

cat(
  sprintf(
    "P = %.4f\n",
    P_cont
  )
)

cat("\nSchoenfeld test:\n")
print(
  ph_cont
)


###############################################################################
# 12. Median-based KSCC score groups
###############################################################################

median_score <- median(
  surv_df$KSCC_Score,
  na.rm = TRUE
)

surv_df <- surv_df %>%
  mutate(
    KSCC_group =
      ifelse(
        KSCC_Score >= median_score,
        "High KSCC score",
        "Low KSCC score"
      ),
    KSCC_group =
      factor(
        KSCC_group,
        levels = c(
          "Low KSCC score",
          "High KSCC score"
        )
      )
  )

cat(
  "\nMedian KSCC score =",
  median_score,
  "\n"
)

cat("\nGroup sizes:\n")
print(
  table(
    surv_df$KSCC_group
  )
)

cat("\nDeaths/censored by group:\n")
print(
  table(
    surv_df$KSCC_group,
    surv_df$OS_status
  )
)


###############################################################################
# 13. Kaplan-Meier and log-rank analysis
###############################################################################

km_fit <- survfit(
  Surv(
    OS_time,
    OS_status
  ) ~ KSCC_group,
  data = surv_df
)

logrank_test <- survdiff(
  Surv(
    OS_time,
    OS_status
  ) ~ KSCC_group,
  data = surv_df
)

logrank_chisq <- logrank_test$chisq

logrank_p <- 1 - pchisq(
  logrank_chisq,
  df = 1
)

cat(
  "\n============================================================\n"
)
cat("LOG-RANK TEST\n")
cat(
  "============================================================\n"
)

cat(
  sprintf(
    "Chi-square = %.3f\n",
    logrank_chisq
  )
)

cat(
  sprintf(
    "P = %.4f\n",
    logrank_p
  )
)


###############################################################################
# 14. High vs low KSCC score Cox model
###############################################################################

cox_group <- coxph(
  Surv(
    OS_time,
    OS_status
  ) ~ KSCC_group,
  data = surv_df
)

cox_group_sum <- summary(
  cox_group
)

HR_group <- cox_group_sum$coefficients[
  1,
  "exp(coef)"
]

P_group <- cox_group_sum$coefficients[
  1,
  "Pr(>|z|)"
]

CI_group <- exp(
  confint(
    cox_group
  )
)

cat(
  "\n============================================================\n"
)
cat("HIGH VS LOW KSCC SCORE\n")
cat(
  "============================================================\n"
)

cat(
  sprintf(
    "High vs Low HR = %.3f\n",
    HR_group
  )
)

cat(
  sprintf(
    "95%% CI = %.3f - %.3f\n",
    CI_group[1],
    CI_group[2]
  )
)

cat(
  sprintf(
    "P = %.4f\n",
    P_group
  )
)


###############################################################################
# 15. Kaplan-Meier figure
###############################################################################

km_plot <- ggsurvplot(
  km_fit,
  data = surv_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Overall survival time (days)",
  ylab = "Overall survival probability",
  legend.title = "",
  legend.labs = c(
    "Low KSCC score",
    "High KSCC score"
  ),
  risk.table.height = 0.25,
  ggtheme = theme_classic(
    base_size = 14
  )
)

pdf(
  "Supplementary_Figure_TCGA_LUSC_KSCC_Score_OS_KM.pdf",
  width = 7,
  height = 7
)

print(
  km_plot
)

dev.off()


###############################################################################
# 16. Save survival results
###############################################################################

survival_results <- data.frame(
  Analysis = c(
    "Continuous KSCC score per 1 SD",
    "High vs Low KSCC score",
    "Log-rank test"
  ),
  HR = c(
    HR_cont,
    HR_group,
    NA_real_
  ),
  CI_lower = c(
    CI_cont[1],
    CI_group[1],
    NA_real_
  ),
  CI_upper = c(
    CI_cont[2],
    CI_group[2],
    NA_real_
  ),
  P_value = c(
    P_cont,
    P_group,
    logrank_p
  ),
  stringsAsFactors = FALSE
)

write.csv(
  survival_results,
  "TCGA_LUSC_KSCC_Score_Survival_Results.csv",
  row.names = FALSE
)


###############################################################################
# 17. Final compact summary
###############################################################################

tp63_result <- correlation_results %>%
  filter(Gene == "TP63")

krt6a_result <- correlation_results %>%
  filter(Gene == "KRT6A")

krt15_result <- correlation_results %>%
  filter(Gene == "KRT15")

cat(
  "\n============================================================\n"
)
cat("FINAL TCGA SUMMARY\n")
cat(
  "============================================================\n"
)

cat(
  "Primary tumor samples:",
  nrow(tcga_score_df),
  "\n"
)

cat(
  sprintf(
    "TP63: rho = %.3f, raw P = %.4g, BH-adjusted P = %.4g\n",
    tp63_result$Spearman_rho,
    tp63_result$P_value,
    tp63_result$Adjusted_P_value
  )
)

cat(
  sprintf(
    "KRT6A: rho = %.3f, raw P = %.4g, BH-adjusted P = %.4g\n",
    krt6a_result$Spearman_rho,
    krt6a_result$P_value,
    krt6a_result$Adjusted_P_value
  )
)

cat(
  sprintf(
    "KRT15: rho = %.3f, raw P = %.4g, BH-adjusted P = %.4g\n",
    krt15_result$Spearman_rho,
    krt15_result$P_value,
    krt15_result$Adjusted_P_value
  )
)

cat(
  sprintf(
    "Survival matched N = %d; deaths = %d\n",
    nrow(surv_df),
    sum(surv_df$OS_status == 1)
  )
)

cat(
  sprintf(
    "Continuous score: HR = %.3f (95%% CI %.3f-%.3f), P = %.4f\n",
    HR_cont,
    CI_cont[1],
    CI_cont[2],
    P_cont
  )
)

cat(
  sprintf(
    "High vs Low: HR = %.3f (95%% CI %.3f-%.3f), P = %.4f\n",
    HR_group,
    CI_group[1],
    CI_group[2],
    P_group
  )
)

cat(
  sprintf(
    "Log-rank: chi-square = %.3f, P = %.4f\n",
    logrank_chisq,
    logrank_p
  )
)


###############################################################################
# Optional:
# sessionInfo()
###############################################################################

# sessionInfo()


###############################################################################
# End of TCGA_analysis.R
###############################################################################
