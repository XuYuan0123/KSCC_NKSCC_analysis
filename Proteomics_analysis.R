###############################################################################
# Proteomics_analysis.R
#
# Proteomic analysis for:
# "Prognostic Implications and Proteomic Classification of
#  Keratinizing versus Non-Keratinizing Lung Squamous Cell Carcinoma"
#
# Required input files:
#   1. data/KSCC_proteomics.xlsx
#   2. data/NKSCC_proteomic.xlsx
#
# Main verified analysis:
# 5697 raw proteins
# -> missingness <20%
# -> 5293 proteins
# -> unified KNN imputation (k = 10)
# -> log2(x + 1)
# -> quantile normalization
# -> limma, eBayes(trend = TRUE)
# -> 410 DE proteins
# -> PLS-DA
# -> mean VIP > 1: 152 proteins
# -> LASSO: 35 proteins
# -> Boruta: 30 proteins
# -> intersection: 10 proteins
#
# TCGA analyses are provided separately in TCGA_analysis.R.
###############################################################################


###############################################################################
# 0. Packages
###############################################################################

library(readxl)
library(tibble)
library(dplyr)
library(impute)
library(preprocessCore)
library(limma)
library(mixOmics)
library(glmnet)
library(Boruta)
library(caret)
library(pROC)


###############################################################################
# 1. Read raw proteomic data
###############################################################################

KSCC <- read_excel(
  "data/KSCC_proteomics.xlsx"
)

NKSCC <- read_excel(
  "data/NKSCC_proteomics.xlsx"
)


# Confirm identical protein identities and order
stopifnot(
  nrow(KSCC) == nrow(NKSCC),
  identical(
    rownames(KSCC),
    rownames(NKSCC)
  )
)


# Combine raw matrices
combined_expr <- cbind(
  KSCC,
  NKSCC
)

storage.mode(
  combined_expr
) <- "numeric"


n_KSCC <- ncol(KSCC)
n_NKSCC <- ncol(NKSCC)


sample_group <- c(
  rep(
    "KSCC",
    n_KSCC
  ),
  rep(
    "NKSCC",
    n_NKSCC
  )
)


# Basic checks
stopifnot(
  nrow(combined_expr) == 5697,
  ncol(combined_expr) == 97,
  n_KSCC == 67,
  n_NKSCC == 30
)


cat(
  "\n============================================================\n"
)

cat(
  "RAW DATA\n"
)

cat(
  "============================================================\n"
)

cat(
  "Proteins:",
  nrow(combined_expr),
  "\n"
)

cat(
  "Samples :",
  ncol(combined_expr),
  "\n"
)

cat(
  "KSCC    :",
  n_KSCC,
  "\n"
)

cat(
  "NKSCC   :",
  n_NKSCC,
  "\n"
)


###############################################################################
# 2. Global preprocessing for the primary discovery analysis
#
# IMPORTANT:
# The verified manuscript-generating branch performs KNN imputation after
# combining all 97 samples. Do not replace this with group-wise KNN.
###############################################################################

full_data <- data.frame(
  SampleID = colnames(
    combined_expr
  ),
  Group = sample_group,
  t(
    combined_expr
  )
)


# Missingness filtering
na_count <- colSums(
  is.na(
    full_data
  )
)


full_data_clean <- full_data[
  ,
  na_count /
    nrow(full_data) <
    0.20,
  drop = FALSE
]


stopifnot(
  nrow(full_data_clean) == 97,
  ncol(full_data_clean) == 5295
)


# Extract expression matrix
expr_matrix <- full_data_clean[
  ,
  -c(1, 2),
  drop = FALSE
] %>%
  as.matrix()


mode(
  expr_matrix
) <- "numeric"


# protein × sample
expr_matrix <- t(
  expr_matrix
)


stopifnot(
  nrow(expr_matrix) == 5293,
  ncol(expr_matrix) == 97
)


cat(
  "\nProteins after missingness filtering:",
  nrow(expr_matrix),
  "\n"
)

cat(
  "Missing values before KNN:",
  sum(
    is.na(
      expr_matrix
    )
  ),
  "\n"
)


###############################################################################
# 3. Unified KNN imputation
###############################################################################

if (
  sum(
    is.na(
      expr_matrix
    )
  ) > 0
) {

  expr_imputed <- impute::impute.knn(
    expr_matrix,
    k = 10
  )$data

} else {

  expr_imputed <- expr_matrix
}


stopifnot(
  sum(
    is.na(
      expr_imputed
    )
  ) == 0
)


cat(
  "Missing values after KNN:",
  sum(
    is.na(
      expr_imputed
    )
  ),
  "\n"
)


###############################################################################
# 4. log2 transformation and quantile normalization
###############################################################################

expr_log2 <- log2(
  expr_imputed + 1
)


expr_normalized <- limma::normalizeQuantiles(
  expr_log2
)


rownames(
  expr_normalized
) <- rownames(
  expr_log2
)

colnames(
  expr_normalized
) <- colnames(
  expr_log2
)


stopifnot(
  nrow(expr_normalized) == 5293,
  ncol(expr_normalized) == 97
)


cat(
  "Normalized matrix:",
  dim(expr_normalized),
  "\n"
)


###############################################################################
# 5. Differential protein analysis: limma
###############################################################################

group_main <- factor(
  full_data_clean$Group,
  levels = c(
    "KSCC",
    "NKSCC"
  )
)


design <- model.matrix(
  ~ 0 + group_main
)


colnames(
  design
) <- levels(
  group_main
)


fit <- limma::lmFit(
  expr_normalized,
  design
)


contrast_matrix <- limma::makeContrasts(
  KSCC_vs_NKSCC =
    KSCC - NKSCC,
  levels = design
)


fit2 <- limma::contrasts.fit(
  fit,
  contrast_matrix
)


fit2 <- limma::eBayes(
  fit2,
  trend = TRUE
)


diff_results <- limma::topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH"
)


kscc_up <- diff_results %>%
  filter(
    adj.P.Val < 0.05,
    logFC > 1
  ) %>%
  arrange(
    desc(logFC)
  )


nkscc_up <- diff_results %>%
  filter(
    adj.P.Val < 0.05,
    logFC < -1
  ) %>%
  arrange(
    logFC
  )


sig_protein_ids <- rownames(
  diff_results[
    diff_results$adj.P.Val < 0.05 &
      abs(
        diff_results$logFC
      ) > 1,
    ,
    drop = FALSE
  ]
)


stopifnot(
  nrow(diff_results) == 5293,
  nrow(kscc_up) == 194,
  nrow(nkscc_up) == 216,
  length(sig_protein_ids) == 410
)


cat(
  "\n============================================================\n"
)

cat(
  "LIMMA RESULTS\n"
)

cat(
  "============================================================\n"
)

cat(
  "Total proteins analysed:",
  nrow(diff_results),
  "\n"
)

cat(
  "KSCC-up proteins:",
  nrow(kscc_up),
  "\n"
)

cat(
  "NKSCC-up proteins:",
  nrow(nkscc_up),
  "\n"
)

cat(
  "Total DE proteins:",
  length(sig_protein_ids),
  "\n"
)


###############################################################################
# 6. PLS-DA and VIP screening
###############################################################################

plsda_input <- t(
  expr_normalized[
    sig_protein_ids,
    ,
    drop = FALSE
  ]
)


group_factor <- factor(
  group_main,
  levels = c(
    "KSCC",
    "NKSCC"
  )
)


set.seed(123)


plsda_res <- mixOmics::plsda(
  X = plsda_input,
  Y = group_factor,
  ncomp = 2
)


vip_scores <- mixOmics::vip(
  plsda_res
)


vip_df <- data.frame(
  Protein = rownames(
    vip_scores
  ),
  VIP_comp1 =
    vip_scores[, 1],
  VIP_comp2 =
    vip_scores[, 2],
  VIP_mean =
    rowMeans(
      vip_scores[
        ,
        1:2,
        drop = FALSE
      ]
    ),
  stringsAsFactors = FALSE
) %>%
  arrange(
    desc(VIP_mean)
  )


vip_ids <- vip_df$Protein[
  vip_df$VIP_mean > 1
]


input_protein_ids <- intersect(
  sig_protein_ids,
  vip_ids
)


stopifnot(
  length(sig_protein_ids) == 410,
  length(vip_ids) == 152,
  length(input_protein_ids) == 152
)


cat(
  "\n============================================================\n"
)

cat(
  "PLS-DA / VIP RESULTS\n"
)

cat(
  "============================================================\n"
)

cat(
  "Input DE proteins:",
  length(sig_protein_ids),
  "\n"
)

cat(
  "Mean VIP > 1:",
  length(vip_ids),
  "\n"
)


###############################################################################
# 7. LASSO feature selection
###############################################################################

X_all <- t(
  expr_normalized[
    input_protein_ids,
    ,
    drop = FALSE
  ]
)


y_all <- factor(
  group_main,
  levels = c(
    "NKSCC",
    "KSCC"
  )
)


set.seed(123)


cv_fit <- glmnet::cv.glmnet(
  x = as.matrix(
    X_all
  ),
  y = y_all,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "deviance"
)


tmp_coef <- coef(
  cv_fit,
  s = "lambda.min"
)


lasso_final_df <- data.frame(
  Protein =
    rownames(
      tmp_coef
    ),
  Coefficient =
    as.numeric(
      tmp_coef
    ),
  stringsAsFactors = FALSE
) %>%
  filter(
    Protein != "(Intercept)",
    Coefficient != 0
  )


lasso_ids <- lasso_final_df$Protein


stopifnot(
  length(lasso_ids) == 35
)


cat(
  "\n============================================================\n"
)

cat(
  "LASSO RESULTS\n"
)

cat(
  "============================================================\n"
)

cat(
  "Input proteins:",
  ncol(X_all),
  "\n"
)

cat(
  "lambda.min:",
  cv_fit$lambda.min,
  "\n"
)

cat(
  "Selected proteins:",
  length(lasso_ids),
  "\n"
)


###############################################################################
# 8. Boruta feature selection
###############################################################################

set.seed(123)


boruta_res <- Boruta::Boruta(
  x = X_all,
  y = y_all,
  pValue = 0.01,
  mcAdj = TRUE,
  maxRuns = 500
)


boruta_ids <- Boruta::getSelectedAttributes(
  boruta_res,
  withTentative = FALSE
)


stopifnot(
  length(boruta_ids) == 30
)


cat(
  "\n============================================================\n"
)

cat(
  "BORUTA RESULTS\n"
)

cat(
  "============================================================\n"
)

cat(
  "Input proteins:",
  ncol(X_all),
  "\n"
)

cat(
  "Confirmed proteins:",
  length(boruta_ids),
  "\n"
)


###############################################################################
# 9. Consensus 10-protein panel
###############################################################################

core_proteins <- intersect(
  lasso_ids,
  boruta_ids
)


expected_core <- c(
  "ERGIC2",
  "PKP1",
  "PRXL2B",
  "IDO1",
  "LGALS7B",
  "SIGLEC1",
  "COL5A1",
  "MRPL2",
  "GAMT",
  "CPNE2"
)


stopifnot(
  length(core_proteins) == 10,
  setequal(
    core_proteins,
    expected_core
  )
)


# Keep manuscript order
final_signature_proteins <-
  expected_core


cat(
  "\n============================================================\n"
)

cat(
  "FINAL FEATURE PANEL\n"
)

cat(
  "============================================================\n"
)

cat(
  "LASSO:",
  length(lasso_ids),
  "\n"
)

cat(
  "Boruta:",
  length(boruta_ids),
  "\n"
)

cat(
  "Intersection:",
  length(core_proteins),
  "\n"
)

print(
  final_signature_proteins
)


###############################################################################
# 10. Fixed 10-protein classifier
###############################################################################

model_x <- as.data.frame(
  t(
    expr_normalized[
      final_signature_proteins,
      ,
      drop = FALSE
    ]
  )
)


rownames(
  model_x
) <- colnames(
  expr_normalized
)


model_data_caret <- model_x


# caret twoClassSummary uses the first factor level as the event class
model_data_caret$Group <- factor(
  as.character(
    group_main
  ),
  levels = c(
    "KSCC",
    "NKSCC"
  )
)


###############################################################################
# 11. Repeated 10-fold CV, 10 repeats
###############################################################################

set.seed(123)


ctrl <- caret::trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction =
    caret::twoClassSummary,
  savePredictions = "final"
)


cv_model <- caret::train(
  Group ~ .,
  data = model_data_caret,
  method = "glm",
  family = binomial,
  metric = "ROC",
  trControl = ctrl
)


# Average repeated out-of-fold probabilities for each sample
cv_pred <- cv_model$pred %>%
  group_by(
    rowIndex,
    obs
  ) %>%
  summarise(
    KSCC =
      mean(KSCC),
    NKSCC =
      mean(NKSCC),
    .groups = "drop"
  ) %>%
  arrange(
    rowIndex
  )


cv_pred$obs <- factor(
  cv_pred$obs,
  levels = c(
    "NKSCC",
    "KSCC"
  )
)


stopifnot(
  nrow(cv_pred) == 97
)


roc_cv <- pROC::roc(
  response =
    cv_pred$obs,
  predictor =
    cv_pred$KSCC,
  levels = c(
    "NKSCC",
    "KSCC"
  ),
  direction = "<",
  quiet = TRUE
)


auc_cv <- as.numeric(
  pROC::auc(
    roc_cv
  )
)


ci_cv <- pROC::ci.auc(
  roc_cv,
  method = "delong"
)


cat(
  "\n============================================================\n"
)

cat(
  "REPEATED 10-FOLD CV × 10\n"
)

cat(
  "============================================================\n"
)

cat(
  "Samples:",
  nrow(cv_pred),
  "\n"
)

cat(
  "AUC:",
  auc_cv,
  "\n"
)

cat(
  "95% CI:",
  as.numeric(ci_cv[1]),
  "-",
  as.numeric(ci_cv[3]),
  "\n"
)


###############################################################################
# 12. Cross-validated KSCC score
###############################################################################

score_summary <- cv_pred %>%
  group_by(
    obs
  ) %>%
  summarise(
    n = n(),
    Mean_KSCC_score =
      mean(KSCC),
    SD_KSCC_score =
      sd(KSCC),
    .groups = "drop"
  )


score_wilcox <- wilcox.test(
  KSCC ~ obs,
  data = cv_pred,
  exact = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "CROSS-VALIDATED KSCC SCORE\n"
)

cat(
  "============================================================\n"
)

print(
  score_summary
)

cat(
  "Wilcoxon P =",
  score_wilcox$p.value,
  "\n"
)


###############################################################################
# 13. Calibration data
###############################################################################

cal_df <- cv_pred %>%
  mutate(
    Event =
      ifelse(
        obs == "KSCC",
        1,
        0
      ),
    PredictedProb =
      KSCC,
    Bin =
      dplyr::ntile(
        PredictedProb,
        5
      )
  ) %>%
  group_by(
    Bin
  ) %>%
  summarise(
    Mean_Predicted =
      mean(PredictedProb),
    Observed =
      mean(Event),
    n = n(),
    .groups = "drop"
  )


cat(
  "\n============================================================\n"
)

cat(
  "CALIBRATION DATA\n"
)

cat(
  "============================================================\n"
)

print(
  cal_df
)


###############################################################################
# 14. Youden cutoff and classification metrics
###############################################################################

best_coords <- pROC::coords(
  roc_cv,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity",
    "accuracy"
  ),
  transpose = FALSE
)


best_cutoff <- as.numeric(
  best_coords$threshold
)


cv_pred$PredictedClass <- ifelse(
  cv_pred$KSCC >=
    best_cutoff,
  "KSCC",
  "NKSCC"
)


cv_pred$PredictedClass <- factor(
  cv_pred$PredictedClass,
  levels = c(
    "NKSCC",
    "KSCC"
  )
)


conf_res <- caret::confusionMatrix(
  data =
    cv_pred$PredictedClass,
  reference =
    cv_pred$obs,
  positive =
    "KSCC"
)


classification_metrics <- data.frame(
  Metric = c(
    "Youden cutoff",
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "Positive predictive value",
    "Negative predictive value",
    "Balanced accuracy",
    "Kappa"
  ),
  Value = c(
    best_cutoff,
    unname(
      conf_res$overall[
        "Accuracy"
      ]
    ),
    unname(
      conf_res$byClass[
        "Sensitivity"
      ]
    ),
    unname(
      conf_res$byClass[
        "Specificity"
      ]
    ),
    unname(
      conf_res$byClass[
        "Pos Pred Value"
      ]
    ),
    unname(
      conf_res$byClass[
        "Neg Pred Value"
      ]
    ),
    unname(
      conf_res$byClass[
        "Balanced Accuracy"
      ]
    ),
    unname(
      conf_res$overall[
        "Kappa"
      ]
    )
  )
)


cat(
  "\n============================================================\n"
)

cat(
  "YOUDEN / CLASSIFICATION RESULTS\n"
)

cat(
  "============================================================\n"
)

print(
  classification_metrics
)


###############################################################################
# 15. Exploratory simpler-model comparison
#
# The exact repeated-CV folds of the 10-protein analysis are reused.
###############################################################################

hist_data <- model_x


hist_data$Group <- factor(
  as.character(
    group_main
  ),
  levels = c(
    "KSCC",
    "NKSCC"
  )
)


hist_index <-
  cv_model$control$index

hist_indexOut <-
  cv_model$control$indexOut


ctrl_benchmark <- caret::trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  index = hist_index,
  indexOut = hist_indexOut,
  classProbs = TRUE,
  summaryFunction =
    caret::twoClassSummary,
  savePredictions = "final"
)


run_hist_benchmark <- function(
    signature,
    model_name
) {

  dat <- hist_data[
    ,
    c(
      signature,
      "Group"
    ),
    drop = FALSE
  ]


  set.seed(123)


  fit_benchmark <- caret::train(
    Group ~ .,
    data = dat,
    method = "glm",
    family = binomial,
    metric = "ROC",
    trControl =
      ctrl_benchmark
  )


  pred_avg <- fit_benchmark$pred %>%
    group_by(
      rowIndex,
      obs
    ) %>%
    summarise(
      KSCC =
        mean(KSCC),
      NKSCC =
        mean(NKSCC),
      .groups = "drop"
    ) %>%
    arrange(
      rowIndex
    )


  pred_avg$obs <- factor(
    pred_avg$obs,
    levels = c(
      "NKSCC",
      "KSCC"
    )
  )


  roc_obj <- pROC::roc(
    response =
      pred_avg$obs,
    predictor =
      pred_avg$KSCC,
    levels = c(
      "NKSCC",
      "KSCC"
    ),
    direction = "<",
    quiet = TRUE
  )


  ci_obj <- pROC::ci.auc(
    roc_obj,
    method = "delong"
  )


  result <- data.frame(
    Model =
      model_name,
    N_proteins =
      length(signature),
    Proteins =
      paste(
        signature,
        collapse = " + "
      ),
    AUC =
      as.numeric(
        pROC::auc(
          roc_obj
        )
      ),
    CI_lower =
      as.numeric(
        ci_obj[1]
      ),
    CI_upper =
      as.numeric(
        ci_obj[3]
      ),
    stringsAsFactors = FALSE
  )


  return(
    list(
      result =
        result,
      fit =
        fit_benchmark,
      predictions =
        pred_avg,
      roc =
        roc_obj
    )
  )
}


benchmark_1 <- run_hist_benchmark(
  "IDO1",
  "1-protein"
)


benchmark_2 <- run_hist_benchmark(
  c(
    "IDO1",
    "PRXL2B"
  ),
  "2-protein"
)


benchmark_3 <- run_hist_benchmark(
  c(
    "IDO1",
    "PRXL2B",
    "SIGLEC1"
  ),
  "3-protein"
)


benchmark_10 <- run_hist_benchmark(
  final_signature_proteins,
  "10-protein"
)


benchmark_results <- bind_rows(
  benchmark_1$result,
  benchmark_2$result,
  benchmark_3$result,
  benchmark_10$result
)


cat(
  "\n============================================================\n"
)

cat(
  "EXPLORATORY MODEL BENCHMARK\n"
)

cat(
  "============================================================\n"
)

print(
  benchmark_results
)


###############################################################################
# 16. Bootstrap OOB internal validation
#
# set.seed(123) is explicitly fixed here for public reproducibility.
###############################################################################

model_data_glm <- model_x


model_data_glm$Group <- factor(
  ifelse(
    group_main == "KSCC",
    "KSCC",
    "NKSCC"
  ),
  levels = c(
    "NKSCC",
    "KSCC"
  )
)


boot_data <- model_data_glm


boot_data$Event <- ifelse(
  boot_data$Group ==
    "KSCC",
  1,
  0
)


feature_cols <- setdiff(
  colnames(
    boot_data
  ),
  c(
    "Group",
    "Event"
  )
)


set.seed(123)


B <- 1000

n <- nrow(
  boot_data
)


boot_auc <- rep(
  NA_real_,
  B
)


for (
  b in seq_len(B)
) {

  boot_idx <- sample(
    seq_len(n),
    size = n,
    replace = TRUE
  )


  oob_idx <- setdiff(
    seq_len(n),
    unique(
      boot_idx
    )
  )


  if (
    length(oob_idx) < 10
  ) {
    next
  }


  if (
    length(
      unique(
        boot_data$Group[
          oob_idx
        ]
      )
    ) < 2
  ) {
    next
  }


  fit_boot <- try(
    glm(
      Event ~ .,
      data =
        boot_data[
          boot_idx,
          c(
            feature_cols,
            "Event"
          ),
          drop = FALSE
        ],
      family =
        binomial
    ),
    silent = TRUE
  )


  if (
    inherits(
      fit_boot,
      "try-error"
    )
  ) {
    next
  }


  pred_oob <- try(
    predict(
      fit_boot,
      newdata =
        boot_data[
          oob_idx,
          feature_cols,
          drop = FALSE
        ],
      type =
        "response"
    ),
    silent = TRUE
  )


  if (
    inherits(
      pred_oob,
      "try-error"
    )
  ) {
    next
  }


  roc_boot <- try(
    pROC::roc(
      response =
        boot_data$Group[
          oob_idx
        ],
      predictor =
        as.numeric(
          pred_oob
        ),
      levels = c(
        "NKSCC",
        "KSCC"
      ),
      direction = "<",
      quiet = TRUE
    ),
    silent = TRUE
  )


  if (
    inherits(
      roc_boot,
      "try-error"
    )
  ) {
    next
  }


  boot_auc[b] <- as.numeric(
    pROC::auc(
      roc_boot
    )
  )
}


boot_auc_clean <- boot_auc[
  !is.na(
    boot_auc
  )
]


boot_auc_mean <- mean(
  boot_auc_clean
)


boot_auc_ci <- quantile(
  boot_auc_clean,
  probs = c(
    0.025,
    0.975
  )
)


bootstrap_auc_summary <- data.frame(
  Mean_AUC =
    boot_auc_mean,
  CI_lower =
    unname(
      boot_auc_ci[1]
    ),
  CI_upper =
    unname(
      boot_auc_ci[2]
    ),
  Valid_iterations =
    length(
      boot_auc_clean
    )
)


cat(
  "\n============================================================\n"
)

cat(
  "BOOTSTRAP OOB INTERNAL VALIDATION\n"
)

cat(
  "============================================================\n"
)

print(
  bootstrap_auc_summary
)


###############################################################################
# 17. Final full-cohort 10-protein logistic model
#
# These coefficients are used in the separate TCGA workflow.
###############################################################################

final_fit_kscc <- glm(
  Group ~ .,
  data =
    model_data_glm,
  family =
    binomial
)


model_data_glm$KSCC_Score <- predict(
  final_fit_kscc,
  type = "response"
)


full_score_summary <- model_data_glm %>%
  group_by(
    Group
  ) %>%
  summarise(
    n = n(),
    Mean_KSCC_score =
      mean(
        KSCC_Score
      ),
    SD_KSCC_score =
      sd(
        KSCC_Score
      ),
    .groups = "drop"
  )


full_score_wilcox <- wilcox.test(
  KSCC_Score ~ Group,
  data =
    model_data_glm,
  exact = FALSE
)


coef_df_kscc <- data.frame(
  Variable =
    names(
      coef(
        final_fit_kscc
      )
    ),
  Coefficient =
    as.numeric(
      coef(
        final_fit_kscc
      )
    ),
  stringsAsFactors = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "FINAL FULL-COHORT MODEL\n"
)

cat(
  "============================================================\n"
)

print(
  full_score_summary
)

cat(
  "Wilcoxon P =",
  full_score_wilcox$p.value,
  "\n"
)

print(
  coef_df_kscc
)


###############################################################################
# 18. Fully nested 10-fold cross-validation
#
# Final reviewer-response implementation.
#
# Within every outer training fold:
#
# raw protein matrix
# -> training-only missingness filtering
# -> KNN imputation
# -> log2(x + 1)
# -> training-derived quantile normalization
# -> limma
# -> PLS-DA / mean VIP > 1
# -> LASSO
# -> Boruta
# -> LASSO ∩ Boruta
# -> logistic regression
# -> outer held-out prediction
#
# Historical finalized nested-CV details retained:
#   eBayes(fit2) in the nested pipeline
#   inner LASSO nfolds = 5
#   no-consensus-feature fold is retained using training prevalence
###############################################################################

expr_raw <- combined_expr


storage.mode(
  expr_raw
) <- "numeric"


group_nested <- factor(
  c(
    rep(
      "KSCC",
      n_KSCC
    ),
    rep(
      "NKSCC",
      n_NKSCC
    )
  ),
  levels = c(
    "NKSCC",
    "KSCC"
  )
)


names(
  group_nested
) <- colnames(
  expr_raw
)


stopifnot(
  nrow(expr_raw) == 5697,
  ncol(expr_raw) == 97,
  identical(
    names(group_nested),
    colnames(expr_raw)
  )
)


###############################################################################
# 18.1 Helper:
# KNN-impute each held-out sample against the corresponding training matrix
###############################################################################

impute_test_from_train <- function(
    train_mat,
    test_mat,
    k = 10
) {

  test_imp <- matrix(
    NA_real_,
    nrow =
      nrow(test_mat),
    ncol =
      ncol(test_mat),
    dimnames =
      dimnames(test_mat)
  )


  for (
    j in seq_len(
      ncol(test_mat)
    )
  ) {

    combined_mat <- cbind(
      train_mat,
      test_mat[
        ,
        j,
        drop = FALSE
      ]
    )


    imp_combined <- impute::impute.knn(
      combined_mat,
      k = k
    )$data


    test_imp[
      ,
      j
    ] <- imp_combined[
      ,
      ncol(
        imp_combined
      )
    ]
  }


  return(
    test_imp
  )
}


###############################################################################
# 18.2 Create stratified outer folds
###############################################################################

set.seed(123)


outer_folds <- caret::createFolds(
  group_nested,
  k = 10,
  list = TRUE,
  returnTrain = FALSE
)


nested_diagnostics <-
  data.frame()

nested_predictions <-
  data.frame()

nested_selected_features <-
  list()


###############################################################################
# 18.3 Outer loop
###############################################################################

for (
  i in seq_along(
    outer_folds
  )
) {

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "NESTED OUTER FOLD:",
    i,
    "\n"
  )


  test_index <-
    outer_folds[[i]]


  train_index <- setdiff(
    seq_len(
      ncol(expr_raw)
    ),
    test_index
  )


  train_expr <- expr_raw[
    ,
    train_index,
    drop = FALSE
  ]


  test_expr <- expr_raw[
    ,
    test_index,
    drop = FALSE
  ]


  train_group <- factor(
    group_nested[
      train_index
    ],
    levels = c(
      "NKSCC",
      "KSCC"
    )
  )


  test_group <- factor(
    group_nested[
      test_index
    ],
    levels = c(
      "NKSCC",
      "KSCC"
    )
  )


  ###########################################################################
  # Missingness filtering — training data only
  ###########################################################################

  miss_rate <- rowMeans(
    is.na(
      train_expr
    )
  )


  keep_protein <- rownames(
    train_expr
  )[
    miss_rate <= 0.20
  ]


  train_expr <- train_expr[
    keep_protein,
    ,
    drop = FALSE
  ]


  test_expr <- test_expr[
    keep_protein,
    ,
    drop = FALSE
  ]


  n_after_missing <-
    length(
      keep_protein
    )


  ###########################################################################
  # KNN imputation
  ###########################################################################

  train_imp <- impute::impute.knn(
    train_expr,
    k = 10
  )$data


  test_imp <- impute_test_from_train(
    train_mat =
      train_expr,
    test_mat =
      test_expr,
    k = 10
  )


  stopifnot(
    sum(
      is.na(
        train_imp
      )
    ) == 0,
    sum(
      is.na(
        test_imp
      )
    ) == 0
  )


  ###########################################################################
  # log2 transformation
  ###########################################################################

  train_log <- log2(
    train_imp + 1
  )


  test_log <- log2(
    test_imp + 1
  )


  ###########################################################################
  # Quantile normalization:
  # target distribution derived only from the training fold
  ###########################################################################

  q_target <-
    preprocessCore::normalize.quantiles.determine.target(
      train_log
    )


  train_norm <-
    preprocessCore::normalize.quantiles.use.target(
      train_log,
      target =
        q_target
    )


  test_norm <-
    preprocessCore::normalize.quantiles.use.target(
      test_log,
      target =
        q_target
    )


  rownames(
    train_norm
  ) <- rownames(
    train_log
  )


  colnames(
    train_norm
  ) <- colnames(
    train_log
  )


  rownames(
    test_norm
  ) <- rownames(
    test_log
  )


  colnames(
    test_norm
  ) <- colnames(
    test_log
  )


  ###########################################################################
  # limma — training data only
  ###########################################################################

  design_nested <- model.matrix(
    ~ 0 + train_group
  )


  colnames(
    design_nested
  ) <- levels(
    train_group
  )


  fit_nested <- limma::lmFit(
    train_norm,
    design_nested
  )


  contrast_nested <- limma::makeContrasts(
    KSCC - NKSCC,
    levels =
      design_nested
  )


  fit2_nested <- limma::contrasts.fit(
    fit_nested,
    contrast_nested
  )


  # Retained from finalized nested-CV record
  fit2_nested <- limma::eBayes(
    fit2_nested
  )


  diff_nested <- limma::topTable(
    fit2_nested,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )


  sig_nested <- rownames(
    diff_nested[
      diff_nested$adj.P.Val < 0.05 &
        abs(
          diff_nested$logFC
        ) > 1,
      ,
      drop = FALSE
    ]
  )


  n_DE <- length(
    sig_nested
  )


  vip_nested <-
    character(0)

  lasso_nested <-
    character(0)

  boruta_nested <-
    character(0)

  final_nested <-
    character(0)


  ###########################################################################
  # PLS-DA / VIP
  ###########################################################################

  if (
    n_DE >= 2
  ) {

    pls_nested_data <- t(
      train_norm[
        sig_nested,
        ,
        drop = FALSE
      ]
    )


    pls_nested <- mixOmics::plsda(
      X =
        pls_nested_data,
      Y =
        train_group,
      ncomp =
        2
    )


    vip_nested_score <- mixOmics::vip(
      pls_nested
    )


    vip_nested_value <- rowMeans(
      vip_nested_score[
        ,
        1:2,
        drop = FALSE
      ],
      na.rm = TRUE
    )


    vip_nested <- names(
      vip_nested_value[
        vip_nested_value >
          1
      ]
    )
  }


  n_VIP <- length(
    vip_nested
  )


  ###########################################################################
  # LASSO
  ###########################################################################

  if (
    n_VIP > 0
  ) {

    X_lasso_nested <- t(
      train_norm[
        vip_nested,
        ,
        drop = FALSE
      ]
    )


    y_lasso_nested <- ifelse(
      train_group ==
        "KSCC",
      1,
      0
    )


    set.seed(
      1000 + i
    )


    lasso_nested_fit <- tryCatch(

      glmnet::cv.glmnet(
        X_lasso_nested,
        y_lasso_nested,
        family =
          "binomial",
        alpha =
          1,
        nfolds =
          5
      ),

      error = function(e) {
        NULL
      }

    )


    if (
      !is.null(
        lasso_nested_fit
      )
    ) {

      coef_lasso_nested <- coef(
        lasso_nested_fit,
        s =
          "lambda.min"
      )


      lasso_nested <- rownames(
        coef_lasso_nested
      )[
        as.numeric(
          coef_lasso_nested
        ) != 0
      ]


      lasso_nested <- setdiff(
        lasso_nested,
        "(Intercept)"
      )
    }


    #########################################################################
    # Boruta
    #########################################################################

    X_boruta_nested <- as.data.frame(
      t(
        train_norm[
          vip_nested,
          ,
          drop = FALSE
        ]
      )
    )


    set.seed(
      2000 + i
    )


    boruta_nested <- tryCatch(

      {

        boruta_nested_fit <- Boruta::Boruta(
          x =
            X_boruta_nested,
          y =
            train_group,
          pValue =
            0.01,
          mcAdj =
            TRUE,
          maxRuns =
            500,
          doTrace =
            0
        )


        Boruta::getSelectedAttributes(
          boruta_nested_fit,
          withTentative =
            FALSE
        )

      },

      error = function(e) {
        character(0)
      }

    )


    final_nested <- intersect(
      lasso_nested,
      boruta_nested
    )
  }


  n_LASSO <- length(
    lasso_nested
  )


  n_Boruta <- length(
    boruta_nested
  )


  n_final <- length(
    final_nested
  )


  nested_selected_features[[paste0(
    "Fold",
    i
  )]] <- final_nested


  ###########################################################################
  # Logistic regression
  #
  # Zero consensus feature:
  # retain fold and use training KSCC prevalence.
  ###########################################################################

  status <-
    "model_fitted"


  if (
    n_final == 0
  ) {

    pred_nested <- rep(
      mean(
        train_group ==
          "KSCC"
      ),
      length(
        test_index
      )
    )


    status <-
      "no_consensus_feature"

  } else {

    train_model_nested <- as.data.frame(
      t(
        train_norm[
          final_nested,
          ,
          drop = FALSE
        ]
      )
    )


    test_model_nested <- as.data.frame(
      t(
        test_norm[
          final_nested,
          ,
          drop = FALSE
        ]
      )
    )


    safe_names <- make.names(
      final_nested,
      unique = TRUE
    )


    colnames(
      train_model_nested
    ) <- safe_names


    colnames(
      test_model_nested
    ) <- safe_names


    train_model_nested$Outcome <- ifelse(
      train_group ==
        "KSCC",
      1,
      0
    )


    nested_model <- tryCatch(

      suppressWarnings(
        glm(
          Outcome ~ .,
          data =
            train_model_nested,
          family =
            binomial()
        )
      ),

      error = function(e) {
        NULL
      }

    )


    if (
      is.null(
        nested_model
      )
    ) {

      pred_nested <- rep(
        mean(
          train_group ==
            "KSCC"
        ),
        length(
          test_index
        )
      )


      status <-
        "model_failure"

    } else {

      pred_nested <- suppressWarnings(
        predict(
          nested_model,
          newdata =
            test_model_nested,
          type =
            "response"
        )
      )


      bad_pred <- !is.finite(
        pred_nested
      )


      if (
        any(
          bad_pred
        )
      ) {

        pred_nested[
          bad_pred
        ] <- mean(
          train_group ==
            "KSCC"
        )
      }
    }
  }


  ###########################################################################
  # Fold-specific AUC
  ###########################################################################

  if (
    length(
      unique(
        pred_nested
      )
    ) == 1 ||
      length(
        unique(
          test_group
        )
      ) < 2
  ) {

    auc_nested_fold <-
      0.5

  } else {

    roc_nested_fold <- pROC::roc(
      response =
        test_group,
      predictor =
        pred_nested,
      levels = c(
        "NKSCC",
        "KSCC"
      ),
      direction =
        "<",
      quiet =
        TRUE
    )


    auc_nested_fold <- as.numeric(
      pROC::auc(
        roc_nested_fold
      )
    )
  }


  ###########################################################################
  # Save diagnostics
  ###########################################################################

  nested_diagnostics <- rbind(
    nested_diagnostics,
    data.frame(
      Fold =
        i,
      N_train =
        length(
          train_index
        ),
      N_test =
        length(
          test_index
        ),
      Proteins_after_missing =
        n_after_missing,
      DE =
        n_DE,
      VIP =
        n_VIP,
      LASSO =
        n_LASSO,
      Boruta =
        n_Boruta,
      Final =
        n_final,
      AUC =
        auc_nested_fold,
      Status =
        status,
      stringsAsFactors =
        FALSE
    )
  )


  ###########################################################################
  # Save outer held-out predictions
  ###########################################################################

  nested_predictions <- rbind(
    nested_predictions,
    data.frame(
      Fold =
        i,
      Sample =
        colnames(
          test_expr
        ),
      True_Group =
        as.character(
          test_group
        ),
      Prediction =
        as.numeric(
          pred_nested
        ),
      stringsAsFactors =
        FALSE
    )
  )


  cat(
    "Train:",
    length(train_index),
    "| Test:",
    length(test_index),
    "| Proteins:",
    n_after_missing,
    "| DE:",
    n_DE,
    "| VIP:",
    n_VIP,
    "| LASSO:",
    n_LASSO,
    "| Boruta:",
    n_Boruta,
    "| Final:",
    n_final,
    "| AUC:",
    round(
      auc_nested_fold,
      4
    ),
    "| Status:",
    status,
    "\n"
  )
}


###############################################################################
# 19. Critical nested-CV checks
###############################################################################

stopifnot(
  nrow(
    nested_diagnostics
  ) == 10
)


stopifnot(
  nrow(
    nested_predictions
  ) ==
    ncol(
      expr_raw
    )
)


stopifnot(
  length(
    unique(
      nested_predictions$Sample
    )
  ) ==
    ncol(
      expr_raw
    )
)


stopifnot(
  all(
    table(
      nested_predictions$Sample
    ) == 1
  )
)


###############################################################################
# 20. Pooled nested out-of-fold ROC/AUC
###############################################################################

nested_predictions$True_Group <- factor(
  nested_predictions$True_Group,
  levels = c(
    "NKSCC",
    "KSCC"
  )
)


roc_nested_oof <- pROC::roc(
  response =
    nested_predictions$True_Group,
  predictor =
    nested_predictions$Prediction,
  levels = c(
    "NKSCC",
    "KSCC"
  ),
  direction =
    "<",
  quiet =
    TRUE
)


nested_pooled_auc <- as.numeric(
  pROC::auc(
    roc_nested_oof
  )
)


nested_pooled_ci <- pROC::ci.auc(
  roc_nested_oof,
  method = "delong"
)


nested_mean_fold_auc <- mean(
  nested_diagnostics$AUC
)


nested_sd_fold_auc <- sd(
  nested_diagnostics$AUC
)


nested_feature_frequency <- sort(
  table(
    unlist(
      nested_selected_features
    )
  ),
  decreasing = TRUE
)


nested_feature_frequency_df <- data.frame(
  Protein =
    names(
      nested_feature_frequency
    ),
  Frequency =
    as.numeric(
      nested_feature_frequency
    ),
  stringsAsFactors =
    FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "FULLY NESTED 10-FOLD CV\n"
)

cat(
  "============================================================\n"
)


print(
  nested_diagnostics
)


cat(
  "\nMean fold-wise AUC:",
  nested_mean_fold_auc,
  "\n"
)


cat(
  "Fold-wise AUC SD:",
  nested_sd_fold_auc,
  "\n"
)


cat(
  "Pooled OOF AUC:",
  nested_pooled_auc,
  "\n"
)


cat(
  "Pooled OOF 95% CI:",
  as.numeric(
    nested_pooled_ci[1]
  ),
  "-",
  as.numeric(
    nested_pooled_ci[3]
  ),
  "\n"
)


cat(
  "Outer folds retained:",
  nrow(
    nested_diagnostics
  ),
  "\n"
)


cat(
  "Held-out samples:",
  nrow(
    nested_predictions
  ),
  "/",
  ncol(
    expr_raw
  ),
  "\n"
)


###############################################################################
# 21. Final compact summary
###############################################################################

final_summary <- data.frame(
  Item = c(
    "Raw proteins",
    "Proteins after missingness filtering",
    "Differential proteins",
    "VIP candidates",
    "LASSO proteins",
    "Boruta proteins",
    "Consensus proteins",
    "Repeated-CV AUC",
    "Repeated-CV CI lower",
    "Repeated-CV CI upper",
    "Youden cutoff",
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "Balanced accuracy",
    "Kappa",
    "Bootstrap mean AUC",
    "Bootstrap CI lower",
    "Bootstrap CI upper",
    "Nested pooled AUC",
    "Nested pooled CI lower",
    "Nested pooled CI upper"
  ),
  Value = c(
    nrow(
      combined_expr
    ),
    nrow(
      expr_normalized
    ),
    length(
      sig_protein_ids
    ),
    length(
      vip_ids
    ),
    length(
      lasso_ids
    ),
    length(
      boruta_ids
    ),
    length(
      core_proteins
    ),
    auc_cv,
    as.numeric(
      ci_cv[1]
    ),
    as.numeric(
      ci_cv[3]
    ),
    best_cutoff,
    unname(
      conf_res$overall[
        "Accuracy"
      ]
    ),
    unname(
      conf_res$byClass[
        "Sensitivity"
      ]
    ),
    unname(
      conf_res$byClass[
        "Specificity"
      ]
    ),
    unname(
      conf_res$byClass[
        "Balanced Accuracy"
      ]
    ),
    unname(
      conf_res$overall[
        "Kappa"
      ]
    ),
    boot_auc_mean,
    unname(
      boot_auc_ci[1]
    ),
    unname(
      boot_auc_ci[2]
    ),
    nested_pooled_auc,
    as.numeric(
      nested_pooled_ci[1]
    ),
    as.numeric(
      nested_pooled_ci[3]
    )
  ),
  stringsAsFactors = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "FINAL SUMMARY\n"
)

cat(
  "============================================================\n"
)


print(
  final_summary,
  row.names = FALSE
)


###############################################################################
# Optional:
# sessionInfo()
###############################################################################

# sessionInfo()


###############################################################################
# End of Proteomics_analysis.R
###############################################################################
