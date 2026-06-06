# Chronic Kidney Disease Risk Factor Analysis using R

## Overview

This project investigates the clinical factors associated with Chronic Kidney Disease (CKD) and develops predictive models using binary and ordinal logistic regression. The objective is to identify which biomarkers : hemoglobin, serum creatinine, blood urea, albumin, and GFR ; most strongly predict both CKD status and disease severity, and to evaluate model performance in a screening context.

This is a follow-up to the [Diabetes Risk Factor Analysis](https://github.com/adharsh-github/diabetes-risk-analysis-r), applying a more advanced modelling pipeline to a richer, discretised clinical dataset.

---

## Dataset

The analysis is based on a discretised variant of the UCI CKD Dataset (`ckd-dataset-v2.csv`), containing 200 patients across 29 clinical variables, including:

- Hemoglobin (g/dL)
- Serum Creatinine (mg/dL)
- Blood Urea (mg/dL)
- Albumin (0–5 ordinal scale)
- Blood Glucose Random (mg/dL)
- GFR (mL/min/1.73m²)
- CKD Stage (s1–s5, ordered severity)
- Binary flags: hypertension, diabetes mellitus, anaemia, pedal oedema, red blood cells, pus cells, etc.
- Target: `class` (ckd / notckd)

**Target distribution:** 128 CKD (64%) · 72 Non-CKD (36%) · 200 complete cases

---

## Methodology

### Data Preprocessing

- Dropped metadata rows; standardised missing value codes (`""`, `"?"`, `"\t?"`)
- Encoded all clinical variables as ordered factors with clinically meaningful level labels
- Binary clinical flags converted to unordered factors

### Exploratory Analysis

- Stage distribution bar chart (s1–s5)
- Hemoglobin distribution by CKD status
- Serum creatinine distribution by CKD status
- Hemoglobin progression across CKD stages (bubble chart, CKD patients only)
- GFR distribution across CKD stages (bubble chart, confirming dataset consistency)

### Statistical Testing

- Chi-square and Fisher's exact tests (simulated p-values, B = 10,000) for association between each clinical variable and CKD class
- Fisher's exact used where expected cell counts fell below 5

### Predictive Modelling

**Binary logistic regression** (CKD vs non-CKD):
- Predictors: hemoglobin, blood urea, albumin (numeric rank encoding)
- Interpreted via odds ratios with 95% confidence intervals

**Ordinal logistic regression** (CKD stage s1 → s5):
- Restricted to CKD patients; modelled severity progression
- Predictors: hemoglobin, serum creatinine, blood urea, albumin
- Fitted using `MASS::polr()` with logistic link

### Model Evaluation

- Confusion matrix at threshold 0.50
- Accuracy, sensitivity, specificity
- ROC curve and AUC (`pROC`)
- Threshold sweep (0.10–0.90) to visualise sensitivity–specificity trade-off

---

## Results

| Metric | Value |
|---|---|
| AUC | **0.989** |
| Accuracy (threshold 0.50) | 95.2% |
| Sensitivity | 96.6% |
| Specificity | 93.1% |

### Chi-Square / Fisher's Test

| Variable | p-value | Significance |
|---|---|---|
| hemo | 0.0001 | *** |
| sc | < 0.0001 | *** |
| bu | 0.0001 | *** |
| al | 0.0001 | *** |
| bgr | 0.0001 | *** |
| htn | < 0.0001 | *** |
| dm | < 0.0001 | *** |
| ane | < 0.0001 | *** |
| pe | < 0.0001 | *** |
| grf | 0.1287 | — |

### Key Findings

- Hemoglobin is the strongest individual classifier of CKD status (OR = 0.113 per rank step, p < 0.001) ; each rank increase in hemoglobin reduces CKD odds by ~89%
- Blood urea nitrogen dominates CKD severity staging (ordinal logistic t ≈ 23.5), far outweighing hemoglobin and albumin
- Albumin exhibits complete separation in the dataset, inflating its binary logistic OR to implausible values ; a known artefact in small perfectly-separating predictors
- GFR was the only variable not significantly associated with CKD class (p = 0.13), likely due to sparse cell counts across 9 GFR bins triggering Fisher's exact test
- The model achieves a substantial improvement over the prior diabetes analysis (AUC 0.989 vs 0.83)

---

## Visualisations

All plots are saved to the `images/` folder:

| File | Description |
|---|---|
| `01_stage_distribution.png` | Patient counts across CKD stages s1–s5 |
| `02_hemoglobin_vs_class.png` | Hemoglobin distribution by CKD status |
| `03_creatinine_vs_class.png` | Serum creatinine distribution by CKD status |
| `04_hemo_across_stages.png` | Hemoglobin progression across CKD stages (bubble chart) |
| `05_grf_across_stages.png` | GFR distribution across CKD stages (bubble chart) |
| `06_roc_curve.png` | ROC curve with AUC annotation |
| `07_threshold_tradeoff.png` | Sensitivity vs specificity vs accuracy across thresholds |

---

## Limitations

- Albumin's complete separation in the binary logistic model produces an uninterpretable OR (~4 × 10⁸); the predictor should be regularised or removed for coefficient interpretation
- Ordinal model drops 41 observations due to missingness; imputation was not applied
- Dataset is relatively small (200 patients); results should be validated on a larger cohort
- GFR's non-significance is likely a dataset artefact (sparse cells) rather than a true clinical finding, since GFR is definitionally linked to CKD staging

---

## Conclusion

This analysis confirms hemoglobin as the dominant classifier of CKD presence, and blood urea nitrogen as the key driver of severity progression from stage s1 to s5 — consistent with the pathophysiology of renal anaemia and declining filtration. The binary logistic model achieves near-perfect discrimination (AUC = 0.989), a significant improvement over the AUC of 0.83 in the prior diabetes project. The ordinal regression framework adds clinical value by modelling disease severity rather than just presence, enabling more actionable risk stratification.

---

## Tools & Libraries

- R
- `ggplot2` — visualisation
- `dplyr` — data manipulation
- `MASS` — ordinal logistic regression (`polr`)
- `pROC` — ROC curve and AUC
- `gridExtra` — plot arrangement

---

## Future Work

- Apply missing data imputation (MICE or median imputation) to retain the full 200-observation dataset in ordinal models
- Address complete separation in albumin via Firth's penalised logistic regression (`logistf`)
- Explore tree-based classifiers (Random Forest, XGBoost) for comparison
- Extend ordinal modelling to include GFR and creatinine as continuous covariates
- Validate on an external CKD cohort

---

## Related Projects

- [Diabetes Risk Factor Analysis](https://github.com/adharsh-github/diabetes-risk-analysis-r) — logistic regression and ROC evaluation on the Pima Indians Diabetes Dataset
