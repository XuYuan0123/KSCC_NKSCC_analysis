# Prognostic Implications and Proteomic Classification of Keratinizing versus Non-Keratinizing Lung Squamous Cell Carcinoma

This repository contains the R code used for the main statistical and bioinformatic analyses in the study:

**“Prognostic Implications and Proteomic Classification of Keratinizing versus Non-Keratinizing Lung Squamous Cell Carcinoma.”**

The analyses include population-based survival analysis, institutional-cohort validation, proteomic feature discovery and classifier development, and TCGA-LUSC transcriptomic biological relevance assessment.

## Repository structure

```text
README.md
SEER_analysis.R
Inhouse_analysis.R
Proteomics_analysis.R
TCGA_analysis.R
```

### `SEER_analysis.R`

Performs the SEER-based prognostic analyses, including:

- cohort preprocessing;
- baseline comparisons;
- Kaplan–Meier survival analysis;
- propensity score matching;
- Cox proportional hazards regression;
- proportional hazards assumption testing;
- matched-cohort analyses.

The final propensity score matching procedure uses 1:1 nearest-neighbor matching without replacement.

### `Inhouse_analysis.R`

Performs the institutional clinical-cohort analyses, including:

- clinical-variable preprocessing;
- baseline comparisons;
- Kaplan–Meier survival analysis;
- reverse Kaplan–Meier follow-up estimation;
- univariable Cox regression;
- multivariable Cox regression;
- proportional hazards assumption testing.

The institutional patient-level raw data are not included in this repository because they contain non-public clinical information.

### `Proteomics_analysis.R`

Performs the proteomic discovery and classifier-development workflow, including:

- proteomic data import and preprocessing;
- missing-value filtering;
- K-nearest-neighbor imputation;
- log2 transformation;
- quantile normalization;
- limma differential-expression analysis;
- partial least-squares discriminant analysis;
- VIP-based candidate selection;
- LASSO feature selection;
- Boruta feature selection;
- derivation of the final 10-protein panel;
- logistic-regression classifier construction;
- repeated 10-fold cross-validation;
- Youden-index-based operating-point assessment;
- bootstrap internal validation;
- fully nested 10-fold cross-validation;
- exploratory comparison with simpler protein models.

The final 10-protein panel is:

```text
ERGIC2
PKP1
PRXL2B
IDO1
LGALS7B
SIGLEC1
COL5A1
MRPL2
GAMT
CPNE2
```

The finalized preprocessing and feature-selection sequence is:

```text
5,697 raw proteins
→ 5,293 proteins after missingness filtering
→ 410 differentially expressed proteins
→ 152 proteins with mean VIP > 1
→ 35 proteins selected by LASSO
→ 30 proteins confirmed by Boruta
→ 10 overlapping core proteins
```

For the final LASSO analysis, `lambda.min = 0.008405748`.

### `TCGA_analysis.R`

Performs the TCGA-LUSC transcriptomic biological relevance and survival analyses, including:

- analysis of 501 TCGA-LUSC primary tumor samples;
- extraction of transcripts corresponding to the 10-protein panel;
- log2(TPM + 1) transformation;
- within-cohort Z-score standardization;
- calculation of a TCGA-derived KSCC score using coefficients from the proteomic logistic model;
- Spearman correlation with `TP63`;
- Spearman correlation with eight keratinization-related keratin genes:
  `KRT5`, `KRT6A`, `KRT6B`, `KRT13`, `KRT14`, `KRT15`, `KRT16`, and `KRT17`;
- Benjamini–Hochberg correction across the nine prespecified correlation tests;
- exploratory overall-survival analysis.

Because RNA-seq expression and LFQ proteomic abundance are measured on different scales, the TCGA-derived KSCC score is interpreted as an **ordinal molecular index rather than a calibrated KSCC probability**. The TCGA analysis is therefore used to assess biological relevance and is not treated as formal external validation of classification performance.

## Software

The analyses were conducted in R. The manuscript reports use of:

```text
R version 4.2.0
```

Major R packages used across the scripts include:

```text
survival
survminer
MatchIt
tableone
limma
impute
preprocessCore
mixOmics
glmnet
Boruta
caret
pROC
dplyr
tidyr
tibble
ggplot2
ggpubr
pheatmap
readxl
```

Additional package dependencies are loaded directly within the individual scripts.

## Input data

The repository contains analysis code only. Raw datasets are not redistributed where access restrictions, patient confidentiality, or source-data agreements apply.

### SEER cohort

SEER patient-level data should be obtained through the appropriate SEER access process. Local file paths in `SEER_analysis.R` should be adapted to the user's own SEER data location.

### Institutional clinical cohort

The institutional patient-level dataset is not publicly distributed because it contains non-public clinical information. Researchers seeking access should follow the data-availability statement and institutional requirements described in the manuscript.

### Proteomic cohort

The proteomic data used in this study are not publicly included in this repository because they are being used in ongoing research.

`Proteomics_analysis.R` contains the complete preprocessing, differential-expression, feature-selection, and classifier-development workflow used in the study. Researchers with access to the corresponding proteomic expression matrices may adapt the input paths in the script to reproduce the analysis.
### TCGA-LUSC

`TCGA_analysis.R` expects the following files:

```text
gencode.v36.annotation.gtf.gene.probemap
TCGA-LUSC.star_tpm.tsv.gz
TCGA-LUSC.survival.tsv
```

The transcriptomic analyses use 501 TCGA-LUSC primary tumor samples. Survival analyses are restricted to patients with available survival information.

## Reproducibility notes

Several stochastic procedures use an explicit random seed to improve reproducibility. Where specified in the scripts, the principal seed is:

```r
set.seed(123)
```

The reported bootstrap internal-validation result corresponds to the reproducible seeded analysis:

```text
Mean AUC = 0.814
95% CI = 0.629–0.958
```

The repeated 10-fold cross-validation result for the fixed 10-protein model is:

```text
AUC = 0.826
95% CI = 0.731–0.921
```

The fully nested 10-fold cross-validation result is:

```text
Pooled AUC = 0.719
95% CI = 0.603–0.835
```

These validation procedures address different questions. Repeated cross-validation evaluates the fixed 10-protein panel, whereas the fully nested analysis repeats preprocessing, feature selection, and model fitting within each outer training fold and therefore provides a more stringent estimate of pipeline-level out-of-sample performance.

## Running the analyses

Each script is designed to run independently after the required input data are available and local file paths have been configured. There is no required execution order.

The TCGA script contains the finalized coefficients from the clean 10-protein proteomic logistic-regression model so that it can be run independently without loading an external RDS model object.

## Data interpretation

The 10-protein model was developed for molecular discrimination between keratinizing and non-keratinizing lung squamous cell carcinoma in the discovery proteomic cohort. Internal validation results should be interpreted in the context of the study design and cohort size.

The TCGA analysis evaluates transcriptomic biological relevance of the proteomics-derived score. It should not be interpreted as direct cross-platform validation of the proteomic classifier or as evidence that the TCGA-derived score is a calibrated subtype probability.

## Citation

If you use this code in academic work, please cite the corresponding article. Citation details will be updated upon publication.
