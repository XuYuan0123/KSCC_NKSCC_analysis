# COMPLETE SEER ANALYSIS
# Run from repository root.


###############################################################################
# 01_SEER_data_preparation.R
#
# Purpose:
#   Recreate the analysis-ready SEER cohort used in the manuscript.
#
# Input:
#   data/SEER_export.csv
#
# Output:
#   results/SEER_clean_analysis_dataset.rds
#   results/SEER_clean_analysis_dataset.csv
###############################################################################

required_packages <- c("dplyr", "stringr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

input_file <- file.path("data", "SEER_export.csv")
if (!file.exists(input_file)) {
  stop(
    "SEER input file not found: ", input_file,
    "\nPlace the authorized SEER case-listing export at data/SEER_export.csv."
  )
}

# check.names = TRUE reproduces the syntactic column names used in the original
# analysis script after read.csv().
seer_raw <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = TRUE
)

required_columns <- c(
  "First.malignant.primary.indicator",
  "Diagnostic.Confirmation",
  "Survival.months.flag",
  "Survival.months",
  "Vital.status.recode..study.cutoff.used.",
  "ICD.O.3.Hist.behav",
  "Race.recode..White..Black..Other.",
  "Grade.Recode..thru.2017.",
  "Grade.Pathological..2018..",
  "Derived.EOD.2018.Stage.Group.Recode..2018..",
  "Derived.AJCC.Stage.Group..7th.ed..2010.2015.",
  "Derived.AJCC.Stage.Group..6th.ed..2004.2015.",
  "Derived.EOD.2018.T.Recode..2018..",
  "Derived.AJCC.T..7th.ed..2010.2015.",
  "Derived.AJCC.T..6th.ed..2004.2015.",
  "Derived.EOD.2018.N.Recode..2018..",
  "Derived.AJCC.N..7th.ed..2010.2015.",
  "Derived.AJCC.N..6th.ed..2004.2015.",
  "Derived.EOD.2018.M.Recode..2018..",
  "Derived.AJCC.M..7th.ed..2010.2015.",
  "Derived.AJCC.M..6th.ed..2004.2015.",
  "RX.Summ..Surg.Prim.Site..1998..",
  "Chemotherapy.recode..yes..no.unk.",
  "Radiation.recode",
  "Age.recode.with.single.ages.and.90.",
  "Marital.status.at.diagnosis",
  "Sex"
)

missing_columns <- setdiff(required_columns, names(seer_raw))
if (length(missing_columns) > 0) {
  stop(
    "The SEER export is missing required columns:\n",
    paste0(" - ", missing_columns, collapse = "\n"),
    "\n\nThe variable labels may differ in your SEER*Stat release/export."
  )
}

# -----------------------------------------------------------------------------
# 1. Analysis-level eligibility filters
# -----------------------------------------------------------------------------

dat <- seer_raw |>
  dplyr::filter(
    First.malignant.primary.indicator == "Yes",
    Diagnostic.Confirmation == "Positive histology",
    Survival.months.flag ==
      "Complete dates are available and there are more than 0 days of survival"
  ) |>
  dplyr::filter(
    as.character(Survival.months) != "0000",
    as.character(Survival.months) != "Unknown"
  ) |>
  dplyr::mutate(
    Survival_Months = suppressWarnings(
      as.numeric(as.character(Survival.months))
    ),
    Status_OS = ifelse(
      Vital.status.recode..study.cutoff.used. == "Dead", 1L, 0L
    )
  ) |>
  dplyr::filter(
    !is.na(Survival_Months),
    Survival_Months > 0
  )

# -----------------------------------------------------------------------------
# 2. Histological subtype
#    KSCC  = ICD-O-3 8071/3
#    NKSCC = ICD-O-3 8072/3 or 8073/3
# -----------------------------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    Subtype = dplyr::case_when(
      grepl("8071", ICD.O.3.Hist.behav) ~ "Keratinizing",
      grepl("8072|8073", ICD.O.3.Hist.behav) ~ "Non-keratinizing",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(Subtype))

# -----------------------------------------------------------------------------
# 3. Race
#    This reproduces the grouping used in the original analysis:
#    White, Black, and all remaining values grouped as Other.
# -----------------------------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    Race_Group = dplyr::case_when(
      Race.recode..White..Black..Other. == "White" ~ "White",
      Race.recode..White..Black..Other. == "Black" ~ "Black",
      TRUE ~ "Other"
    )
  )

# -----------------------------------------------------------------------------
# 4. Histological grade
#    Prefer the 2018+ pathological-grade variable; otherwise use the
#    pre-2018 grade recode.
# -----------------------------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    Grade_Final = dplyr::case_when(
      Grade.Pathological..2018.. == "1" ~ "G1",
      Grade.Pathological..2018.. == "2" ~ "G2",
      Grade.Pathological..2018.. == "3" ~ "G3",
      Grade.Pathological..2018.. == "4" ~ "G4",
      Grade.Recode..thru.2017. == "Well differentiated; Grade I" ~ "G1",
      Grade.Recode..thru.2017. == "Moderately differentiated; Grade II" ~ "G2",
      Grade.Recode..thru.2017. == "Poorly differentiated; Grade III" ~ "G3",
      Grade.Recode..thru.2017. ==
        "Undifferentiated; anaplastic; Grade IV" ~ "G4",
      TRUE ~ "Unknown"
    )
  )

# -----------------------------------------------------------------------------
# 5. Overall TNM stage
#    Priority: EOD 2018 stage -> AJCC 7th -> AJCC 6th.
#    Summary Stage is not used to back-fill overall TNM stage.
# -----------------------------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    Stage_Final = dplyr::case_when(
      grepl("^1", Derived.EOD.2018.Stage.Group.Recode..2018..) ~ "I",
      grepl("^2", Derived.EOD.2018.Stage.Group.Recode..2018..) ~ "II",
      grepl("^3", Derived.EOD.2018.Stage.Group.Recode..2018..) ~ "III",
      grepl("^4", Derived.EOD.2018.Stage.Group.Recode..2018..) ~ "IV",

      grepl("^IA|^IB", Derived.AJCC.Stage.Group..7th.ed..2010.2015.) ~ "I",
      grepl("^IIA|^IIB|^II$",
            Derived.AJCC.Stage.Group..7th.ed..2010.2015.) ~ "II",
      grepl("^IIIA|^IIIB|^IIIC|^III$",
            Derived.AJCC.Stage.Group..7th.ed..2010.2015.) ~ "III",
      grepl("^IV", Derived.AJCC.Stage.Group..7th.ed..2010.2015.) ~ "IV",

      grepl("^IA|^IB", Derived.AJCC.Stage.Group..6th.ed..2004.2015.) ~ "I",
      grepl("^IIA|^IIB",
            Derived.AJCC.Stage.Group..6th.ed..2004.2015.) ~ "II",
      grepl("^IIIA|^IIIB",
            Derived.AJCC.Stage.Group..6th.ed..2004.2015.) ~ "III",
      grepl("^IV", Derived.AJCC.Stage.Group..6th.ed..2004.2015.) ~ "IV",

      TRUE ~ "Unknown"
    )
  )

# -----------------------------------------------------------------------------
# 6. T stage
# -----------------------------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    T_Final = dplyr::case_when(
      grepl("^T1", Derived.EOD.2018.T.Recode..2018..) ~ "T1",
      grepl("^T2", Derived.EOD.2018.T.Recode..2018..) ~ "T2",
      grepl("^T3", Derived.EOD.2018.T.Recode..2018..) ~ "T3",
      grepl("^T4", Derived.EOD.2018.T.Recode..2018..) ~ "T4",

      grepl("^T1", Derived.AJCC.T..7th.ed..2010.2015.) ~ "T1",
      grepl("^T2", Derived.AJCC.T..7th.ed..2010.2015.) ~ "T2",
      grepl("^T3", Derived.AJCC.T..7th.ed..2010.2015.) ~ "T3",
      grepl("^T4", Derived.AJCC.T..7th.ed..2010.2015.) ~ "T4",

      grepl("^T1", Derived.AJCC.T..6th.ed..2004.2015.) ~ "T1",
      grepl("^T2", Derived.AJCC.T..6th.ed..2004.2015.) ~ "T2",
      grepl("^T3", Derived.AJCC.T..6th.ed..2004.2015.) ~ "T3",
      grepl("^T4", Derived.AJCC.T..6th.ed..2004.2015.) ~ "T4",

      TRUE ~ "Unknown"
    )
  )

# -----------------------------------------------------------------------------
# 7. N stage
# -----------------------------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    N_Final = dplyr::case_when(
      Derived.EOD.2018.N.Recode..2018.. %in%
        c("N0", "N1", "N2", "N3") ~
        as.character(Derived.EOD.2018.N.Recode..2018..),

      Derived.AJCC.N..7th.ed..2010.2015. %in%
        c("N0", "N1", "N2", "N3") ~
        as.character(Derived.AJCC.N..7th.ed..2010.2015.),

      Derived.AJCC.N..6th.ed..2004.2015. %in%
        c("N0", "N1", "N2", "N3") ~
        as.character(Derived.AJCC.N..6th.ed..2004.2015.),

      TRUE ~ "Unknown"
    )
  )

# -----------------------------------------------------------------------------
# 8. M stage
# -----------------------------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    M_Final = dplyr::case_when(
      Derived.EOD.2018.M.Recode..2018.. == "M0" ~ "M0",
      grepl("^M1", Derived.EOD.2018.M.Recode..2018..) ~ "M1",

      Derived.AJCC.M..7th.ed..2010.2015. == "M0" ~ "M0",
      grepl("^M1", Derived.AJCC.M..7th.ed..2010.2015.) ~ "M1",

      Derived.AJCC.M..6th.ed..2004.2015. == "M0" ~ "M0",
      Derived.AJCC.M..6th.ed..2004.2015. == "M1" ~ "M1",

      TRUE ~ "Unknown"
    )
  )

# -----------------------------------------------------------------------------
# 9. Treatment variables
# -----------------------------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    Surgery = dplyr::case_when(
      RX.Summ..Surg.Prim.Site..1998.. == "0" ~ "No/Unknown",
      RX.Summ..Surg.Prim.Site..1998.. %in% c("Blank(s)", "99") ~
        "No/Unknown",
      TRUE ~ "Yes"
    ),

    Chemo = dplyr::case_when(
      Chemotherapy.recode..yes..no.unk. == "Yes" ~ "Yes",
      TRUE ~ "No/Unknown"
    ),

    Radiation = dplyr::case_when(
      Radiation.recode %in% c(
        "Beam radiation",
        "Combination of beam with implants or isotopes",
        "Radioactive implants (includes brachytherapy) (1988+)",
        "Radiation, NOS  method or source not specified"
      ) ~ "Yes",
      Radiation.recode == "Refused (1988+)" ~ "No/Unknown",
      TRUE ~ "No/Unknown"
    )
  )

# -----------------------------------------------------------------------------
# 10. Age and marital status
# -----------------------------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    Age_Numeric = suppressWarnings(
      as.numeric(
        stringr::str_extract(
          Age.recode.with.single.ages.and.90.,
          "\\d+"
        )
      )
    ),
    Age_Group_65 = ifelse(Age_Numeric < 65, "<65", ">=65"),

    Marital_Status = dplyr::case_when(
      Marital.status.at.diagnosis ==
        "Married (including common law)" ~ "Married",
      TRUE ~ "No/Unknown"
    )
  )

# -----------------------------------------------------------------------------
# 11. Complete-key-variable analysis cohort
# -----------------------------------------------------------------------------

data_strict <- dat |>
  dplyr::filter(
    T_Final != "Unknown",
    N_Final != "Unknown",
    M_Final != "Unknown",
    Grade_Final != "Unknown",
    Stage_Final != "Unknown"
  )

# -----------------------------------------------------------------------------
# 12. Explicit factor coding/reference levels
#
# Cox-reference categories match Supplementary Table S1:
# KSCC, age <65, female, Black race, married, grade I, stage I,
# T1, N0, M0, no/unknown treatment.
# -----------------------------------------------------------------------------

data_strict <- data_strict |>
  dplyr::mutate(
    Subtype = factor(
      Subtype,
      levels = c("Keratinizing", "Non-keratinizing")
    ),
    Age_Group_65 = factor(
      Age_Group_65,
      levels = c("<65", ">=65")
    ),
    Sex = factor(
      Sex,
      levels = c("Female", "Male")
    ),
    Race_Group = factor(
      Race_Group,
      levels = c("Black", "White", "Other")
    ),
    Marital_Status = factor(
      Marital_Status,
      levels = c("Married", "No/Unknown")
    ),
    Grade_Final = factor(
      Grade_Final,
      levels = c("G1", "G2", "G3", "G4")
    ),
    Stage_Final = factor(
      Stage_Final,
      levels = c("I", "II", "III", "IV")
    ),
    T_Final = factor(
      T_Final,
      levels = c("T1", "T2", "T3", "T4")
    ),
    N_Final = factor(
      N_Final,
      levels = c("N0", "N1", "N2", "N3")
    ),
    M_Final = factor(
      M_Final,
      levels = c("M0", "M1")
    ),
    Chemo = factor(
      Chemo,
      levels = c("No/Unknown", "Yes")
    ),
    Radiation = factor(
      Radiation,
      levels = c("No/Unknown", "Yes")
    ),
    Surgery = factor(
      Surgery,
      levels = c("No/Unknown", "Yes")
    ),
    Group_ID = ifelse(Subtype == "Keratinizing", 1L, 0L)
  )

# -----------------------------------------------------------------------------
# 13. Manuscript-level verification
# -----------------------------------------------------------------------------

cat("\n===== CLEAN SEER COHORT =====\n")
cat("Total N =", nrow(data_strict), "\n")
print(table(data_strict$Subtype, useNA = "ifany"))
cat("Deaths =", sum(data_strict$Status_OS == 1, na.rm = TRUE), "\n")
cat("Censored =", sum(data_strict$Status_OS == 0, na.rm = TRUE), "\n")

expected_n <- 6397L
expected_kscc <- 3921L
expected_nkscc <- 2476L
expected_deaths <- 3925L

observed_counts <- table(data_strict$Subtype)

if (nrow(data_strict) != expected_n ||
    unname(observed_counts["Keratinizing"]) != expected_kscc ||
    unname(observed_counts["Non-keratinizing"]) != expected_nkscc ||
    sum(data_strict$Status_OS == 1, na.rm = TRUE) != expected_deaths) {
  warning(
    "The cleaned cohort does not match the manuscript counts. ",
    "Check the SEER release/export and case-selection settings."
  )
}

saveRDS(
  data_strict,
  file.path("results", "SEER_clean_analysis_dataset.rds")
)

write.csv(
  data_strict,
  file.path("results", "SEER_clean_analysis_dataset.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


###############################################################################
# 02_SEER_PSM_and_balance.R
#
# Purpose:
#   Reproduce the final main propensity-score model, 1:1 matching,
#   standardized balance assessment, and manuscript Table 1.
#
# Input:
#   results/SEER_clean_analysis_dataset.rds
#
# Output:
#   results/PSM_clinical_model.rds
#   results/SEER_PSM_dataset.rds
#   results/SEER_PSM_dataset.csv
#   results/Table1_SEER_pre_post_PSM.csv
#   results/Table1_SMD_summary.csv
#   results/PSM_balance_cobalt_standardized.csv
###############################################################################

required_packages <- c("dplyr", "MatchIt", "cobalt", "tableone")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

input_file <- file.path("results", "SEER_clean_analysis_dataset.rds")
if (!file.exists(input_file)) {
  stop("Run R/01_SEER_data_preparation.R first.")
}

data_strict <- readRDS(input_file)

# -----------------------------------------------------------------------------
# 1. Main propensity-score matching
#
# Saved final model specification verified from the original MatchIt object:
# Group_ID ~ Age_Group_65 + Sex + Race_Group + Marital_Status +
#            Stage_Final + Surgery + Chemo + Radiation
#
# Important:
# - Grade and T/N/M are NOT separate PSM covariates.
# - Logistic regression estimates propensity scores.
# - 1:1 greedy nearest-neighbor matching without replacement.
# - Caliper = 0.05 SD of the estimated propensity score.
# - m.order = "largest" is explicitly stated; this is the MatchIt default
#   when a propensity score is estimated.
# -----------------------------------------------------------------------------

PSM_clinical_model <- MatchIt::matchit(
  Group_ID ~
    Age_Group_65 +
    Sex +
    Race_Group +
    Marital_Status +
    Stage_Final +
    Surgery +
    Chemo +
    Radiation,
  data = data_strict,
  method = "nearest",
  distance = "glm",
  link = "logit",
  caliper = 0.05,
  std.caliper = TRUE,
  ratio = 1,
  replace = FALSE,
  discard = "none",
  m.order = "largest"
)

SEER_PSM_dataset <- MatchIt::match_data(
  PSM_clinical_model,
  data = data_strict,
  drop.unmatched = TRUE
)

# Ensure Cox reference levels remain exactly as specified in script 01.
SEER_PSM_dataset$Subtype <- factor(
  as.character(SEER_PSM_dataset$Subtype),
  levels = c("Keratinizing", "Non-keratinizing")
)

# -----------------------------------------------------------------------------
# 2. Match verification
# -----------------------------------------------------------------------------

cat("\n===== FINAL PSM MODEL =====\n")
print(PSM_clinical_model$call)

cat("\n===== MATCHED SAMPLE SIZE =====\n")
print(table(SEER_PSM_dataset$Subtype, useNA = "ifany"))
cat("Total matched N =", nrow(SEER_PSM_dataset), "\n")
cat(
  "Matched pairs =",
  length(unique(SEER_PSM_dataset$subclass)),
  "\n"
)

cat("\n===== CALIPER CHECK =====\n")
cat("Stored caliper on propensity-score scale =",
    PSM_clinical_model$caliper, "\n")
cat("SD of estimated propensity score =",
    sd(PSM_clinical_model$distance), "\n")
cat(
  "Stored caliper / SD(PS) =",
  unname(PSM_clinical_model$caliper /
           sd(PSM_clinical_model$distance)),
  "\n"
)

matched_counts <- table(SEER_PSM_dataset$Subtype)
if (nrow(SEER_PSM_dataset) != 4952L ||
    unname(matched_counts["Keratinizing"]) != 2476L ||
    unname(matched_counts["Non-keratinizing"]) != 2476L ||
    length(unique(SEER_PSM_dataset$subclass)) != 2476L) {
  warning(
    "Matched cohort does not reproduce 2,476 matched pairs. ",
    "Check data order, software version, and source SEER export."
  )
}

# -----------------------------------------------------------------------------
# 3. Cobalt standardized balance
# -----------------------------------------------------------------------------

bal_std <- cobalt::bal.tab(
  PSM_clinical_model,
  un = TRUE,
  binary = "std",
  thresholds = c(m = 0.1)
)

balance_detail <- as.data.frame(bal_std$Balance)
balance_detail$Term <- rownames(balance_detail)
rownames(balance_detail) <- NULL
balance_detail <- balance_detail |>
  dplyr::select(Term, dplyr::everything())

write.csv(
  balance_detail,
  file.path("results", "PSM_balance_cobalt_standardized.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -----------------------------------------------------------------------------
# 4. Reconstruct manuscript Table 1
# -----------------------------------------------------------------------------

matching_vars <- c(
  "Age_Group_65",
  "Sex",
  "Race_Group",
  "Marital_Status",
  "Stage_Final",
  "Surgery",
  "Chemo",
  "Radiation"
)

table_vars <- c(
  "Age_Group_65",
  "Sex",
  "Race_Group",
  "Marital_Status",
  "Grade_Final",
  "Stage_Final",
  "T_Final",
  "N_Final",
  "M_Final",
  "Chemo",
  "Radiation",
  "Surgery"
)

display_names <- c(
  Age_Group_65 = "Age",
  Sex = "Sex",
  Race_Group = "Race",
  Marital_Status = "Marital status",
  Grade_Final = "Grade",
  Stage_Final = "TNM stage",
  T_Final = "T stage",
  N_Final = "N stage",
  M_Final = "M stage",
  Chemo = "Chemotherapy",
  Radiation = "Radiation",
  Surgery = "Surgery"
)

display_levels <- list(
  Age_Group_65 = c("<65", ">=65"),
  Sex = c("Male", "Female"),
  Race_Group = c("White", "Black", "Other"),
  Marital_Status = c("Married", "No/Unknown"),
  Grade_Final = c("G1", "G2", "G3", "G4"),
  Stage_Final = c("I", "II", "III", "IV"),
  T_Final = c("T1", "T2", "T3", "T4"),
  N_Final = c("N0", "N1", "N2", "N3"),
  M_Final = c("M0", "M1"),
  Chemo = c("Yes", "No/Unknown"),
  Radiation = c("Yes", "No/Unknown"),
  Surgery = c("Yes", "No/Unknown")
)

format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

get_p_test <- function(x, group) {
  keep <- !is.na(x) & !is.na(group)
  x <- droplevels(as.factor(x[keep]))
  group <- droplevels(as.factor(group[keep]))
  tab <- table(x, group)

  chi <- suppressWarnings(chisq.test(tab))

  if (any(chi$expected < 5)) {
    out <- tryCatch(
      list(p = fisher.test(tab)$p.value, method = "Fisher"),
      error = function(e) {
        set.seed(12345)
        list(
          p = fisher.test(
            tab,
            simulate.p.value = TRUE,
            B = 100000
          )$p.value,
          method = "Fisher"
        )
      }
    )
  } else {
    out <- list(p = chi$p.value, method = "Chi-square")
  }

  out
}

# One overall SMD per categorical factor.
extract_overall_smd <- function(dat, vars) {
  tab1 <- tableone::CreateTableOne(
    vars = vars,
    strata = "Subtype",
    data = dat,
    factorVars = vars,
    test = FALSE
  )
  smd <- tableone::ExtractSmd(tab1)
  if (is.matrix(smd)) smd <- smd[, 1]
  smd <- abs(as.numeric(smd))
  names(smd) <- vars
  smd
}

smd_pre <- extract_overall_smd(data_strict, matching_vars)
smd_post <- extract_overall_smd(SEER_PSM_dataset, matching_vars)

SMD_summary <- data.frame(
  Variable = matching_vars,
  Factor = unname(display_names[matching_vars]),
  SMD_Pre = round(smd_pre[matching_vars], 3),
  SMD_Post = round(smd_post[matching_vars], 3),
  row.names = NULL
)

write.csv(
  SMD_summary,
  file.path("results", "Table1_SMD_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

make_descriptive_table <- function(dat, vars) {
  out <- list()

  # Table 1 displays KSCC first, then NKSCC.
  group <- factor(
    as.character(dat$Subtype),
    levels = c("Keratinizing", "Non-keratinizing"),
    labels = c("KSCC", "NKSCC")
  )

  for (v in vars) {
    levs <- display_levels[[v]]
    x <- factor(as.character(dat[[v]]), levels = levs)

    keep <- !is.na(x) & !is.na(group)
    x2 <- droplevels(x[keep])
    g2 <- droplevels(group[keep])

    tab <- table(x2, g2)
    pct <- prop.table(tab, margin = 2) * 100
    ptest <- get_p_test(x2, g2)

    out[[v]] <- data.frame(
      Variable = v,
      Factor = c(
        unname(display_names[v]),
        rep("", nrow(tab) - 1)
      ),
      Level = rownames(tab),
      KSCC = sprintf(
        "%d (%.1f%%)",
        tab[, "KSCC"],
        pct[, "KSCC"]
      ),
      NKSCC = sprintf(
        "%d (%.1f%%)",
        tab[, "NKSCC"],
        pct[, "NKSCC"]
      ),
      P_value = c(
        format_p(ptest$p),
        rep("", nrow(tab) - 1)
      ),
      Test = c(
        ptest$method,
        rep("", nrow(tab) - 1)
      ),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(out)
}

pre_table <- make_descriptive_table(data_strict, table_vars)
post_table <- make_descriptive_table(SEER_PSM_dataset, table_vars)

Table1_final <- data.frame(
  Variable = pre_table$Variable,
  Factor = pre_table$Factor,
  Level = pre_table$Level,
  Pre_KSCC = pre_table$KSCC,
  Pre_NKSCC = pre_table$NKSCC,
  Pre_P = pre_table$P_value,
  Post_KSCC = post_table$KSCC,
  Post_NKSCC = post_table$NKSCC,
  Post_P = post_table$P_value,
  SMD_Post = "",
  stringsAsFactors = FALSE
)

# The manuscript reports SMD only for variables actually entered in PSM.
# Grade and individual T/N/M are shown descriptively without SMDs.
for (v in table_vars) {
  first_row <- which(Table1_final$Variable == v)[1]

  if (v %in% matching_vars) {
    Table1_final$SMD_Post[first_row] <- sprintf(
      "%.3f",
      smd_post[v]
    )
  } else {
    Table1_final$SMD_Post[first_row] <- "\u2014"
  }
}

write.csv(
  Table1_final,
  file.path("results", "Table1_SEER_pre_post_PSM.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n===== POST-PSM OVERALL SMD =====\n")
print(SMD_summary, row.names = FALSE)

stage_smd <- unname(smd_post["Stage_Final"])
cat("\nPost-PSM overall TNM-stage SMD =", stage_smd, "\n")
if (is.finite(stage_smd) && abs(stage_smd - 0.175) > 0.01) {
  warning(
    "Post-PSM TNM-stage SMD differs materially from the manuscript value (~0.175)."
  )
}

# -----------------------------------------------------------------------------
# 5. Save model/data objects
# -----------------------------------------------------------------------------

saveRDS(
  PSM_clinical_model,
  file.path("results", "PSM_clinical_model.rds")
)

saveRDS(
  SEER_PSM_dataset,
  file.path("results", "SEER_PSM_dataset.rds")
)

write.csv(
  SEER_PSM_dataset,
  file.path("results", "SEER_PSM_dataset.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


###############################################################################
# 03_SEER_survival_analysis.R
#
# Purpose:
#   Reproduce SEER survival analyses and Supplementary Table S1.
#
# Includes:
#   - Reverse Kaplan-Meier follow-up
#   - Kaplan-Meier and log-rank tests before/after PSM
#   - Cluster-robust univariable Cox models in matched pairs
#   - Cluster-robust multivariable Cox model in matched pairs
#   - Schoenfeld proportional-hazards tests
#
# Inputs:
#   results/SEER_clean_analysis_dataset.rds
#   results/SEER_PSM_dataset.rds
###############################################################################

required_packages <- c("survival", "broom", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

pre_file <- file.path("results", "SEER_clean_analysis_dataset.rds")
post_file <- file.path("results", "SEER_PSM_dataset.rds")

if (!file.exists(pre_file) || !file.exists(post_file)) {
  stop("Run scripts 01 and 02 before the survival-analysis script.")
}

data_strict <- readRDS(pre_file)
SEER_PSM_dataset <- readRDS(post_file)

# -----------------------------------------------------------------------------
# 1. Explicit Cox reference levels
# -----------------------------------------------------------------------------

set_survival_levels <- function(dat) {
  dat |>
    dplyr::mutate(
      Subtype = factor(
        as.character(Subtype),
        levels = c("Keratinizing", "Non-keratinizing")
      ),
      Age_Group_65 = factor(
        as.character(Age_Group_65),
        levels = c("<65", ">=65")
      ),
      Sex = factor(
        as.character(Sex),
        levels = c("Female", "Male")
      ),
      Race_Group = factor(
        as.character(Race_Group),
        levels = c("Black", "White", "Other")
      ),
      Marital_Status = factor(
        as.character(Marital_Status),
        levels = c("Married", "No/Unknown")
      ),
      Grade_Final = factor(
        as.character(Grade_Final),
        levels = c("G1", "G2", "G3", "G4")
      ),
      Stage_Final = factor(
        as.character(Stage_Final),
        levels = c("I", "II", "III", "IV")
      ),
      T_Final = factor(
        as.character(T_Final),
        levels = c("T1", "T2", "T3", "T4")
      ),
      N_Final = factor(
        as.character(N_Final),
        levels = c("N0", "N1", "N2", "N3")
      ),
      M_Final = factor(
        as.character(M_Final),
        levels = c("M0", "M1")
      ),
      Chemo = factor(
        as.character(Chemo),
        levels = c("No/Unknown", "Yes")
      ),
      Radiation = factor(
        as.character(Radiation),
        levels = c("No/Unknown", "Yes")
      ),
      Surgery = factor(
        as.character(Surgery),
        levels = c("No/Unknown", "Yes")
      )
    )
}

data_strict <- set_survival_levels(data_strict)
SEER_PSM_dataset <- set_survival_levels(SEER_PSM_dataset)

# -----------------------------------------------------------------------------
# 2. Reverse Kaplan-Meier median follow-up
# -----------------------------------------------------------------------------

reverse_km_all <- survival::survfit(
  survival::Surv(
    Survival_Months,
    1 - Status_OS
  ) ~ 1,
  data = data_strict
)

followup_table <- summary(reverse_km_all)$table
median_followup <- unname(followup_table["median"])
followup_lcl <- unname(followup_table["0.95LCL"])
followup_ucl <- unname(followup_table["0.95UCL"])

cat("\n===== REVERSE KM FOLLOW-UP =====\n")
cat(
  "Median follow-up =",
  median_followup,
  "months (95% CI ",
  followup_lcl, "-",
  followup_ucl, ")\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# 3. Kaplan-Meier and log-rank tests
# -----------------------------------------------------------------------------

km_pre <- survival::survfit(
  survival::Surv(Survival_Months, Status_OS) ~ Subtype,
  data = data_strict
)

logrank_pre <- survival::survdiff(
  survival::Surv(Survival_Months, Status_OS) ~ Subtype,
  data = data_strict
)

km_post <- survival::survfit(
  survival::Surv(Survival_Months, Status_OS) ~ Subtype,
  data = SEER_PSM_dataset
)

logrank_post <- survival::survdiff(
  survival::Surv(Survival_Months, Status_OS) ~ Subtype,
  data = SEER_PSM_dataset
)

logrank_p <- function(x) {
  stats::pchisq(
    x$chisq,
    df = length(x$n) - 1L,
    lower.tail = FALSE
  )
}

km_pre_table <- as.data.frame(summary(km_pre)$table)
km_pre_table$Stratum <- rownames(km_pre_table)
rownames(km_pre_table) <- NULL

km_post_table <- as.data.frame(summary(km_post)$table)
km_post_table$Stratum <- rownames(km_post_table)
rownames(km_post_table) <- NULL

write.csv(
  km_pre_table,
  file.path("results", "SEER_KM_prePSM_summary.csv"),
  row.names = FALSE
)

write.csv(
  km_post_table,
  file.path("results", "SEER_KM_postPSM_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Matched-cohort univariable Cox regression
#
# Every matched-cohort Cox model uses cluster(subclass) to account for
# within-pair dependence. Overall TNM stage is included here (univariable)
# but is not entered into the multivariable model.
# -----------------------------------------------------------------------------

univariable_vars <- c(
  "Subtype",
  "Age_Group_65",
  "Sex",
  "Race_Group",
  "Marital_Status",
  "Grade_Final",
  "Stage_Final",
  "T_Final",
  "N_Final",
  "M_Final",
  "Chemo",
  "Radiation",
  "Surgery"
)

fit_univariable <- function(v) {
  f <- stats::as.formula(
    paste0(
      "survival::Surv(Survival_Months, Status_OS) ~ ",
      v,
      " + cluster(subclass)"
    )
  )

  fit <- survival::coxph(
    f,
    data = SEER_PSM_dataset,
    ties = "efron"
  )

  tt <- broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  )

  tt$Variable <- v
  tt
}

univariable_cox <- dplyr::bind_rows(
  lapply(univariable_vars, fit_univariable)
)

write.csv(
  univariable_cox,
  file.path("results", "SEER_univariable_Cox_cluster_robust.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Matched-cohort multivariable Cox regression
#
# Current manuscript model:
# subtype + age + sex + marital status + grade + T + N + M +
# surgery + chemotherapy + radiotherapy
#
# Race and overall TNM stage are not included in this multivariable model.
# -----------------------------------------------------------------------------

cox_post_adjusted <- survival::coxph(
  survival::Surv(
    Survival_Months,
    Status_OS
  ) ~
    Subtype +
    Age_Group_65 +
    Sex +
    Marital_Status +
    Grade_Final +
    T_Final +
    N_Final +
    M_Final +
    Surgery +
    Chemo +
    Radiation +
    cluster(subclass),
  data = SEER_PSM_dataset,
  ties = "efron"
)

multivariable_cox <- broom::tidy(
  cox_post_adjusted,
  exponentiate = TRUE,
  conf.int = TRUE
)

write.csv(
  multivariable_cox,
  file.path("results", "SEER_multivariable_Cox_cluster_robust.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. Proportional-hazards assumption
# -----------------------------------------------------------------------------

ph_test <- survival::cox.zph(cox_post_adjusted)

ph_table <- as.data.frame(ph_test$table)
ph_table$Term <- rownames(ph_table)
rownames(ph_table) <- NULL
ph_table <- ph_table |>
  dplyr::select(Term, dplyr::everything())

write.csv(
  ph_table,
  file.path("results", "SEER_PH_test.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7. Manuscript-level validation summary
# -----------------------------------------------------------------------------

get_median_by_stratum <- function(km_object, pattern) {
  tab <- summary(km_object)$table
  idx <- grep(pattern, rownames(tab), fixed = TRUE)
  if (length(idx) != 1L) return(NA_real_)
  unname(tab[idx, "median"])
}

subtype_multiv <- multivariable_cox |>
  dplyr::filter(grepl("^Subtype", term))

subtype_univ <- univariable_cox |>
  dplyr::filter(
    Variable == "Subtype",
    grepl("^Subtype", term)
  )

ph_subtype <- ph_table |>
  dplyr::filter(Term == "Subtype")

ph_global <- ph_table |>
  dplyr::filter(Term == "GLOBAL")

validation_summary <- data.frame(
  Metric = c(
    "Original cohort N",
    "Deaths",
    "Matched cohort N",
    "Matched pairs",
    "Median follow-up (months)",
    "Pre-PSM KSCC median OS",
    "Pre-PSM NKSCC median OS",
    "Post-PSM KSCC median OS",
    "Post-PSM NKSCC median OS",
    "Pre-PSM log-rank P",
    "Post-PSM log-rank P",
    "Matched univariable subtype HR",
    "Matched univariable subtype CI low",
    "Matched univariable subtype CI high",
    "Matched multivariable subtype HR",
    "Matched multivariable subtype CI low",
    "Matched multivariable subtype CI high",
    "Matched multivariable subtype P",
    "Schoenfeld subtype P",
    "Schoenfeld GLOBAL P"
  ),
  Value = c(
    nrow(data_strict),
    sum(data_strict$Status_OS == 1),
    nrow(SEER_PSM_dataset),
    length(unique(SEER_PSM_dataset$subclass)),
    median_followup,
    get_median_by_stratum(km_pre, "Subtype=Keratinizing"),
    get_median_by_stratum(km_pre, "Subtype=Non-keratinizing"),
    get_median_by_stratum(km_post, "Subtype=Keratinizing"),
    get_median_by_stratum(km_post, "Subtype=Non-keratinizing"),
    logrank_p(logrank_pre),
    logrank_p(logrank_post),
    subtype_univ$estimate[1],
    subtype_univ$conf.low[1],
    subtype_univ$conf.high[1],
    subtype_multiv$estimate[1],
    subtype_multiv$conf.low[1],
    subtype_multiv$conf.high[1],
    subtype_multiv$p.value[1],
    ph_subtype$p[1],
    ph_global$p[1]
  )
)

write.csv(
  validation_summary,
  file.path("results", "SEER_validation_summary.csv"),
  row.names = FALSE
)

cat("\n===== KEY SEER RESULTS =====\n")
print(validation_summary, row.names = FALSE)

# Soft checks against the manuscript.
check_close <- function(actual, expected, tolerance, label) {
  if (length(actual) == 0L || !is.finite(actual) ||
      abs(actual - expected) > tolerance) {
    warning(
      label,
      " differs from the manuscript check (expected approximately ",
      expected, ")."
    )
  }
}

check_close(median_followup, 51, 0.01, "Median follow-up")
check_close(
  subtype_univ$estimate[1],
  0.88,
  0.02,
  "Matched univariable subtype HR"
)
check_close(
  subtype_multiv$estimate[1],
  0.8014,
  0.01,
  "Matched multivariable subtype HR"
)
check_close(
  ph_subtype$p[1],
  0.9412,
  0.01,
  "Schoenfeld subtype P"
)


###############################################################################
# 04_SEER_figures.R
#
# Purpose:
#   Regenerate the SEER panels corresponding to Figure 2:
#   A. pre-PSM Kaplan-Meier curve
#   B. post-PSM Kaplan-Meier curve
#   C. matched multivariable Cox forest plot
#
# Inputs:
#   results/SEER_clean_analysis_dataset.rds
#   results/SEER_PSM_dataset.rds
###############################################################################

required_packages <- c(
  "survival", "survminer", "ggplot2", "broom", "dplyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

data_strict <- readRDS(
  file.path("results", "SEER_clean_analysis_dataset.rds")
)
SEER_PSM_dataset <- readRDS(
  file.path("results", "SEER_PSM_dataset.rds")
)

# Ensure KSCC is the reference/display-first group.
data_strict$Subtype <- factor(
  as.character(data_strict$Subtype),
  levels = c("Keratinizing", "Non-keratinizing")
)
SEER_PSM_dataset$Subtype <- factor(
  as.character(SEER_PSM_dataset$Subtype),
  levels = c("Keratinizing", "Non-keratinizing")
)

# -----------------------------------------------------------------------------
# Figure 2A: before PSM
# -----------------------------------------------------------------------------

fit_pre <- survival::survfit(
  survival::Surv(Survival_Months, Status_OS) ~ Subtype,
  data = data_strict
)

p_pre <- survminer::ggsurvplot(
  fit_pre,
  data = data_strict,
  xlim = c(0, 96),
  break.time.by = 12,
  conf.int = TRUE,
  conf.int.alpha = 0.15,
  surv.median.line = "hv",
  pval = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  censor.size = 1.5,
  censor.shape = "+",
  size = 1.2,
  palette = c("#D62728", "#1F77B4"),
  legend.labs = c("KSCC", "NKSCC"),
  legend.title = "",
  xlab = "Survival Time (Months)",
  ylab = "Overall Survival Probability",
  ggtheme = ggplot2::theme_classic(base_size = 11)
)

grDevices::pdf(
  file.path("figures", "Figure2A_SEER_KM_prePSM.pdf"),
  width = 8,
  height = 7
)
print(p_pre)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# Figure 2B: after PSM
# -----------------------------------------------------------------------------

fit_post <- survival::survfit(
  survival::Surv(Survival_Months, Status_OS) ~ Subtype,
  data = SEER_PSM_dataset
)

p_post <- survminer::ggsurvplot(
  fit_post,
  data = SEER_PSM_dataset,
  xlim = c(0, 100),
  break.time.by = 20,
  conf.int = TRUE,
  conf.int.alpha = 0.10,
  surv.median.line = "hv",
  pval = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  censor = FALSE,
  size = 1.2,
  palette = c("#EFC000FF", "#0073C2FF"),
  legend.labs = c("KSCC", "NKSCC"),
  legend.title = "",
  xlab = "Survival Time (Months)",
  ylab = "Overall Survival Probability",
  ggtheme = ggplot2::theme_bw(base_size = 11)
)

grDevices::pdf(
  file.path("figures", "Figure2B_SEER_KM_postPSM.pdf"),
  width = 8,
  height = 7
)
print(p_post)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# Figure 2C: matched multivariable Cox forest plot
# -----------------------------------------------------------------------------

# Re-establish reference levels explicitly.
SEER_PSM_dataset <- SEER_PSM_dataset |>
  dplyr::mutate(
    Subtype = factor(Subtype,
                     levels = c("Keratinizing", "Non-keratinizing")),
    Age_Group_65 = factor(Age_Group_65,
                          levels = c("<65", ">=65")),
    Sex = factor(Sex, levels = c("Female", "Male")),
    Marital_Status = factor(Marital_Status,
                            levels = c("Married", "No/Unknown")),
    Grade_Final = factor(Grade_Final,
                         levels = c("G1", "G2", "G3", "G4")),
    T_Final = factor(T_Final,
                     levels = c("T1", "T2", "T3", "T4")),
    N_Final = factor(N_Final,
                     levels = c("N0", "N1", "N2", "N3")),
    M_Final = factor(M_Final,
                     levels = c("M0", "M1")),
    Surgery = factor(Surgery,
                     levels = c("No/Unknown", "Yes")),
    Chemo = factor(Chemo,
                   levels = c("No/Unknown", "Yes")),
    Radiation = factor(Radiation,
                       levels = c("No/Unknown", "Yes"))
  )

cox_fit <- survival::coxph(
  survival::Surv(Survival_Months, Status_OS) ~
    Subtype +
    Age_Group_65 +
    Sex +
    Marital_Status +
    Grade_Final +
    T_Final +
    N_Final +
    M_Final +
    Surgery +
    Chemo +
    Radiation +
    cluster(subclass),
  data = SEER_PSM_dataset,
  ties = "efron"
)

forest_dat <- broom::tidy(
  cox_fit,
  exponentiate = TRUE,
  conf.int = TRUE
)

term_labels <- c(
  "SubtypeNon-keratinizing" = "NKSCC vs KSCC",
  "Age_Group_65>=65" = "Age >=65 vs <65",
  "SexMale" = "Male vs Female",
  "Marital_StatusNo/Unknown" = "Marital: No/Unknown vs Married",
  "Grade_FinalG2" = "Grade II vs I",
  "Grade_FinalG3" = "Grade III vs I",
  "Grade_FinalG4" = "Grade IV vs I",
  "T_FinalT2" = "T2 vs T1",
  "T_FinalT3" = "T3 vs T1",
  "T_FinalT4" = "T4 vs T1",
  "N_FinalN1" = "N1 vs N0",
  "N_FinalN2" = "N2 vs N0",
  "N_FinalN3" = "N3 vs N0",
  "M_FinalM1" = "M1 vs M0",
  "SurgeryYes" = "Surgery: Yes vs No/Unknown",
  "ChemoYes" = "Chemotherapy: Yes vs No/Unknown",
  "RadiationYes" = "Radiotherapy: Yes vs No/Unknown"
)

forest_dat$label <- unname(term_labels[forest_dat$term])
forest_dat$label[is.na(forest_dat$label)] <-
  forest_dat$term[is.na(forest_dat$label)]

forest_dat$label <- factor(
  forest_dat$label,
  levels = rev(forest_dat$label)
)

p_forest <- ggplot2::ggplot(
  forest_dat,
  ggplot2::aes(
    x = estimate,
    y = label
  )
) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    height = 0.15
  ) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    x = "Hazard Ratio (95% CI)",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 11)

ggplot2::ggsave(
  filename = file.path(
    "figures",
    "Figure2C_SEER_multivariable_forest.pdf"
  ),
  plot = p_forest,
  width = 8,
  height = 8
)



writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
