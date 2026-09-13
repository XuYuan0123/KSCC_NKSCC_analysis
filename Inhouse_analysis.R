# ============================================================
# In-house clinical cohort analysis
# Keratinizing vs non-keratinizing lung squamous cell carcinoma
#
# This script reproduces the institutional-cohort survival analyses
# reported in the manuscript.
#
# IMPORTANT:
# - Patient-level institutional data are not included in a public
#   repository because of privacy/institutional restrictions.
# - To reproduce locally, place the authorized de-identified Excel
#   file at the path below, or change input_file accordingly.
# - The original Excel coding is retained:
#       Surgery:      0 = No, 1 = Yes
#       Chemotherapy: 0 = No, 1 = Yes
#       Radiotherapy: 0 = No, 1 = Yes
#       Survival status: 0 = censored/alive, 1 = death
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

required_packages <- c(
  "readxl",
  "survival",
  "survminer",
  "ggplot2",
  "dplyr",
  "broom"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Please install them before running this script."
  )
}

library(readxl)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(broom)


# ------------------------------------------------------------
# 2. Input
# ------------------------------------------------------------

input_file <- "data/inhouse_clinical.xlsx"

# For the original local working file, this can instead be:
# input_file <- "肺鳞癌院内数据.xlsx"

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nPlace the authorized de-identified institutional dataset at this path ",
    "or edit input_file."
  )
}

data_raw <- read_excel(input_file)


# ------------------------------------------------------------
# 3. Required columns and variable construction
# ------------------------------------------------------------

required_columns <- c(
  "组织类型",
  "年龄分组",
  "性别",
  "临床分期",
  "是否进行手术治疗（是=1；否=0）",
  "是否进行化疗（1=是；0=否）",
  "是否进行放疗（是=1；否=0）",
  "生存状态（0=生存；1=死亡）"
)

missing_columns <- setdiff(required_columns, names(data_raw))

if (length(missing_columns) > 0) {
  stop(
    "The following required column(s) are missing:\n",
    paste(missing_columns, collapse = "\n")
  )
}

# The published analysis was verified using the variable "month".
# If "month" is absent, reconstruct it from survival days using 30.44 days/month.
if ("month" %in% names(data_raw)) {
  data_raw$month <- as.numeric(data_raw$month)
} else if ("生存时间（Day）" %in% names(data_raw)) {
  data_raw$month <- as.numeric(data_raw$`生存时间（Day）`) / 30.44
} else {
  stop(
    'No survival-time variable found. The dataset must contain either "month" ',
    'or "生存时间（Day）".'
  )
}

data_raw$`生存状态（0=生存；1=死亡）` <- as.numeric(
  data_raw$`生存状态（0=生存；1=死亡）`
)

# Standardize the few textual variants that occurred in the original file.
stage_chr <- trimws(as.character(data_raw$临床分期))
stage_chr[stage_chr == "Ⅲ"] <- "III"
stage_chr[stage_chr == "Ⅳ"] <- "IV"

age_chr <- trimws(as.character(data_raw$年龄分组))
age_chr[age_chr %in% c("<65", "＜65", "＜65岁")] <- "<65岁"
age_chr[age_chr %in% c(">=65", "≥65", ">=65岁")] <- "≥65岁"

data_inhouse <- data_raw %>%
  mutate(
    Subtype = factor(
      as.character(组织类型),
      levels = c("角化", "非角化"),
      labels = c("KSCC", "NKSCC")
    ),

    Age = factor(
      age_chr,
      levels = c("<65岁", "≥65岁"),
      labels = c("<65", ">=65")
    ),

    Gender = factor(
      as.character(性别),
      levels = c("男", "女"),
      labels = c("Male", "Female")
    ),

    Stage = factor(
      stage_chr,
      levels = c("I", "II", "III", "IV")
    ),

    Surgery = factor(
      as.numeric(`是否进行手术治疗（是=1；否=0）`),
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),

    Chemotherapy = factor(
      as.numeric(`是否进行化疗（1=是；0=否）`),
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),

    Radiotherapy = factor(
      as.numeric(`是否进行放疗（是=1；否=0）`),
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),

    status = as.numeric(`生存状态（0=生存；1=死亡）`)
  )


# ------------------------------------------------------------
# 4. Cohort checks
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("IN-HOUSE COHORT: BASIC INFORMATION\n")
cat("=============================================\n")

cat("Total N =", nrow(data_inhouse), "\n")
cat("Deaths =", sum(data_inhouse$status == 1, na.rm = TRUE), "\n")
cat("Censored/alive =", sum(data_inhouse$status == 0, na.rm = TRUE), "\n")

cat("\nSubtype:\n")
print(table(data_inhouse$Subtype, useNA = "ifany"))

cat("\nAge:\n")
print(table(data_inhouse$Age, useNA = "ifany"))

cat("\nGender:\n")
print(table(data_inhouse$Gender, useNA = "ifany"))

cat("\nTNM stage:\n")
print(table(data_inhouse$Stage, useNA = "ifany"))

cat("\nSurgery (0 = No, 1 = Yes in the original file):\n")
print(table(data_inhouse$Subtype, data_inhouse$Surgery, useNA = "ifany"))

cat("\nChemotherapy (0 = No, 1 = Yes in the original file):\n")
print(table(data_inhouse$Subtype, data_inhouse$Chemotherapy, useNA = "ifany"))

cat("\nRadiotherapy (0 = No, 1 = Yes in the original file):\n")
print(table(data_inhouse$Subtype, data_inhouse$Radiotherapy, useNA = "ifany"))


# Internal validation anchors for the dataset used in the manuscript.
# These checks intentionally stop execution if a different dataset/version
# is supplied, preventing silent mismatch with the published analysis.

stopifnot(
  nrow(data_inhouse) == 283,
  sum(data_inhouse$status == 1, na.rm = TRUE) == 181,
  sum(data_inhouse$status == 0, na.rm = TRUE) == 102,
  unname(table(data_inhouse$Subtype)["KSCC"]) == 65,
  unname(table(data_inhouse$Subtype)["NKSCC"]) == 218
)


# ------------------------------------------------------------
# 5. Baseline characteristics by histological subtype
# ------------------------------------------------------------
# Final manuscript testing rule:
#   Age            -> Pearson's chi-square test
#   Gender         -> Fisher's exact test
#   TNM stage      -> Pearson's chi-square test among patients with
#                     available staging information only
#   Chemotherapy   -> Pearson's chi-square test
#   Radiotherapy   -> Fisher's exact test
#   Surgery        -> Pearson's chi-square test
#
# Fisher's exact P values are reported without a chi-square statistic.

baseline_test_plan <- c(
  Age = "pearson",
  Gender = "fisher",
  Stage = "pearson",
  Chemotherapy = "pearson",
  Radiotherapy = "fisher",
  Surgery = "pearson"
)

cat("\n=============================================\n")
cat("BASELINE CHARACTERISTICS BY SUBTYPE\n")
cat("=============================================\n")

baseline_results <- vector("list", length(baseline_test_plan))
names(baseline_results) <- names(baseline_test_plan)

for (v in names(baseline_test_plan)) {

  cat("\n---", v, "---\n")

  # Complete cases are used for the variable being compared.
  # For Stage, this implements the manuscript analysis among the 221
  # patients with available TNM-stage information.
  keep <- complete.cases(data_inhouse[, c("Subtype", v)])

  tab <- table(
    data_inhouse$Subtype[keep],
    data_inhouse[[v]][keep]
  )

  print(tab)

  if (baseline_test_plan[[v]] == "fisher") {
    test <- fisher.test(tab)
    cat(
      "Fisher's exact test; P =",
      sprintf("%.3f", test$p.value),
      "\n"
    )

    baseline_results[[v]] <- data.frame(
      Variable = v,
      Test = "Fisher's exact test",
      Statistic = NA_real_,
      P_value = unname(test$p.value),
      stringsAsFactors = FALSE
    )

  } else {
    test <- suppressWarnings(chisq.test(tab, correct = FALSE))
    cat(
      "Pearson chi-square =",
      sprintf("%.3f", unname(test$statistic)),
      "; P =",
      sprintf("%.3f", test$p.value),
      "\n"
    )

    baseline_results[[v]] <- data.frame(
      Variable = v,
      Test = "Pearson's chi-square test",
      Statistic = unname(test$statistic),
      P_value = unname(test$p.value),
      stringsAsFactors = FALSE
    )
  }
}

baseline_results <- bind_rows(baseline_results)

cat("\nManuscript-ready baseline test summary:\n")
print(
  baseline_results %>%
    mutate(
      Statistic_report = ifelse(
        is.na(Statistic),
        "—",
        sprintf("%.3f", Statistic)
      ),
      P_report = sprintf("%.3f", P_value)
    ) %>%
    select(Variable, Test, Statistic_report, P_report),
  n = Inf
)

# Reproducibility anchors for Table 2.
# These are checks only; the statistics above are calculated from the raw data.
expected_baseline <- data.frame(
  Variable = c(
    "Age", "Gender", "Stage",
    "Chemotherapy", "Radiotherapy", "Surgery"
  ),
  Statistic = c(3.7638, NA, 5.9187, 2.9015, NA, 5.1303),
  P_value = c(0.05237, 0.7779, 0.1156, 0.08850, 0.1423, 0.02351)
)

for (i in seq_len(nrow(expected_baseline))) {
  v <- expected_baseline$Variable[i]
  got <- baseline_results[baseline_results$Variable == v, ]

  if (!is.na(expected_baseline$Statistic[i])) {
    stopifnot(abs(got$Statistic - expected_baseline$Statistic[i]) < 0.01)
  }
  stopifnot(abs(got$P_value - expected_baseline$P_value[i]) < 0.01)
}

# Additional count checks for the corrected treatment coding in Table 2.
stopifnot(
  unname(table(data_inhouse$Subtype, data_inhouse$Chemotherapy)["KSCC", "Yes"]) == 28,
  unname(table(data_inhouse$Subtype, data_inhouse$Chemotherapy)["NKSCC", "Yes"]) == 69,
  unname(table(data_inhouse$Subtype, data_inhouse$Surgery)["KSCC", "Yes"]) == 13,
  unname(table(data_inhouse$Subtype, data_inhouse$Surgery)["NKSCC", "Yes"]) == 76,
  sum(!is.na(data_inhouse$Stage)) == 221
)


# ------------------------------------------------------------
# 6. Kaplan-Meier overall survival
# ------------------------------------------------------------

fit_km <- survfit(
  Surv(month, status) ~ Subtype,
  data = data_inhouse
)

median_os <- surv_median(fit_km)

cat("\n=============================================\n")
cat("KAPLAN-MEIER OVERALL SURVIVAL\n")
cat("=============================================\n")

print(median_os)

logrank <- survdiff(
  Surv(month, status) ~ Subtype,
  data = data_inhouse
)

logrank_p <- pchisq(
  logrank$chisq,
  df = length(logrank$n) - 1,
  lower.tail = FALSE
)

cat("Log-rank P =", format.pval(logrank_p, digits = 4), "\n")

p_label <- ifelse(
  logrank_p < 0.001,
  "Log-rank P < 0.001",
  paste0("Log-rank P = ", signif(logrank_p, 3))
)

km_plot <- ggsurvplot(
  fit_km,
  data = data_inhouse,
  xlim = c(0, 96),
  break.time.by = 12,
  palette = c("#D62728", "#1F77B4"),
  linetype = "solid",
  size = 1.2,
  censor.size = 1.5,
  censor.shape = "+",
  conf.int = TRUE,
  conf.int.alpha = 0.15,
  surv.median.line = "hv",
  pval = p_label,
  pval.size = 4,
  pval.coord = c(1, 0.15),
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.col = "strata",
  ggtheme = theme_classic(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.2, 0.2),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 9)
    ),
  legend.labs = c("KSCC", "NKSCC"),
  legend.title = "",
  xlab = "Survival Time (Months)",
  ylab = "Overall Survival Probability",
  title = "Institutional cohort"
)

pdf(
  "Figure3_Institutional_Cohort_KM.pdf",
  width = 7.5,
  height = 6.5,
  onefile = FALSE
)
print(km_plot)
dev.off()


# ------------------------------------------------------------
# 7. Reverse Kaplan-Meier median follow-up
# ------------------------------------------------------------

fit_followup <- survfit(
  Surv(month, 1 - status) ~ 1,
  data = data_inhouse
)

followup_table <- summary(fit_followup)$table

median_followup <- as.numeric(followup_table["median"])
followup_lower <- as.numeric(followup_table["0.95LCL"])
followup_upper <- as.numeric(followup_table["0.95UCL"])

cat("\n=============================================\n")
cat("REVERSE KAPLAN-MEIER FOLLOW-UP\n")
cat("=============================================\n")

cat(
  "Median follow-up =",
  round(median_followup, 1),
  "months\n"
)

cat(
  "95% CI =",
  round(followup_lower, 1),
  "to",
  round(followup_upper, 1),
  "months\n"
)


# ------------------------------------------------------------
# 8. Univariable Cox regression
# ------------------------------------------------------------

univariable_variables <- c(
  "Subtype",
  "Gender",
  "Age",
  "Stage",
  "Radiotherapy",
  "Chemotherapy",
  "Surgery"
)

fit_univariable <- function(variable, data) {

  form <- as.formula(
    paste0("Surv(month, status) ~ ", variable)
  )

  fit <- coxph(
    form,
    data = data,
    ties = "efron"
  )

  out <- broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  )

  out$Variable <- variable
  out$Model_N <- fit$n
  out$Events <- fit$nevent

  out
}

univariable_cox <- bind_rows(
  lapply(
    univariable_variables,
    fit_univariable,
    data = data_inhouse
  )
) %>%
  select(
    Variable,
    term,
    Model_N,
    Events,
    estimate,
    conf.low,
    conf.high,
    p.value
  )

cat("\n=============================================\n")
cat("UNIVARIABLE COX REGRESSION\n")
cat("=============================================\n")

print(univariable_cox, n = Inf)


# ------------------------------------------------------------
# 9. Final multivariable Cox regression
# ------------------------------------------------------------
# Variable-selection rule used in the manuscript:
# variables with P < 0.05 in univariable analysis entered the final model;
# histological subtype was forced into the model regardless of univariable P.
#
# Final model:
#   Subtype + Age + overall TNM Stage + Surgery
#
# Complete-case analysis; no imputation.

cox_complete <- data_inhouse %>%
  filter(
    complete.cases(
      month,
      status,
      Subtype,
      Age,
      Stage,
      Surgery
    )
  )

cat("\n=============================================\n")
cat("FINAL MULTIVARIABLE COX MODEL\n")
cat("=============================================\n")

cat("Complete-case N =", nrow(cox_complete), "\n")
cat("Deaths =", sum(cox_complete$status == 1), "\n")

multivariable_cox_model <- coxph(
  Surv(month, status) ~
    Subtype +
    Age +
    Stage +
    Surgery,
  data = cox_complete,
  ties = "efron"
)

print(summary(multivariable_cox_model))

multivariable_cox <- broom::tidy(
  multivariable_cox_model,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  select(
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  )

cat("\nFormatted multivariable estimates:\n")
print(multivariable_cox, n = Inf)


# ------------------------------------------------------------
# 10. Events per regression parameter
# ------------------------------------------------------------

n_events_cox <- multivariable_cox_model$nevent
n_parameters <- length(coef(multivariable_cox_model))
epv <- n_events_cox / n_parameters

cat("\n=============================================\n")
cat("EVENTS PER REGRESSION PARAMETER\n")
cat("=============================================\n")

cat("Events =", n_events_cox, "\n")
cat("Regression parameters =", n_parameters, "\n")
cat("Events per parameter =", round(epv, 2), "\n")


# ------------------------------------------------------------
# 11. Proportional-hazards assumption
# ------------------------------------------------------------

ph_test <- cox.zph(
  multivariable_cox_model,
  transform = "km",
  terms = TRUE
)

cat("\n=============================================\n")
cat("SCHOENFELD RESIDUAL TEST\n")
cat("=============================================\n")

print(ph_test)

ph_table <- as.data.frame(ph_test$table)
ph_table$Variable <- rownames(ph_table)
rownames(ph_table) <- NULL
ph_table <- ph_table[, c("Variable", setdiff(names(ph_table), "Variable"))]


# ------------------------------------------------------------
# 12. Reproducibility checks against the verified manuscript model
# ------------------------------------------------------------

subtype_hr <- unname(
  exp(coef(multivariable_cox_model)["SubtypeNKSCC"])
)

subtype_ci <- exp(
  confint(multivariable_cox_model)["SubtypeNKSCC", ]
)

subtype_p <- summary(multivariable_cox_model)$coefficients[
  "SubtypeNKSCC",
  "Pr(>|z|)"
]

stopifnot(
  nrow(cox_complete) == 221,
  n_events_cox == 132,
  n_parameters == 6,
  abs(epv - 22) < 1e-8,
  abs(subtype_hr - 0.7985) < 0.002,
  abs(subtype_ci[1] - 0.5359) < 0.002,
  abs(subtype_ci[2] - 1.1900) < 0.002,
  abs(subtype_p - 0.26895) < 0.002
)

cat("\n=============================================\n")
cat("VERIFIED ANALYSIS SUMMARY\n")
cat("=============================================\n")

cat("Full cohort N: 283\n")
cat("KSCC: 65; NKSCC: 218\n")
cat("Deaths: 181; censored/alive: 102\n")
cat(
  "Median follow-up:",
  round(median_followup, 1),
  "months (95% CI",
  round(followup_lower, 1),
  "-",
  round(followup_upper, 1),
  ")\n"
)
cat(
  "Median OS: KSCC =",
  round(median_os$median[median_os$strata == "Subtype=KSCC"], 1),
  "months; NKSCC =",
  round(median_os$median[median_os$strata == "Subtype=NKSCC"], 1),
  "months\n"
)
cat("Log-rank P =", format.pval(logrank_p, digits = 4), "\n")
cat("Final Cox complete-case N: 221; deaths: 132\n")
cat(
  "NKSCC vs KSCC adjusted HR =",
  round(subtype_hr, 2),
  "(95% CI",
  round(subtype_ci[1], 2),
  "-",
  round(subtype_ci[2], 2),
  "), P =",
  sprintf("%.3f", subtype_p),
  "\n"
)
cat(
  "Global PH-test P =",
  sprintf("%.3f", ph_test$table["GLOBAL", "p"]),
  "\n"
)
cat("Events per regression parameter =", round(epv, 2), "\n")

cat("\nAll prespecified verification checks passed.\n")
