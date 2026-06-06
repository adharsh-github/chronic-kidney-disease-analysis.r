# =============================================================================
# Chronic Kidney Disease (CKD) Risk Factor Analysis
# Dataset: ckd-dataset-v2.csv (discretised UCI CKD variant)
# Author: Adharsh
# =============================================================================

# ── 1. LOAD LIBRARIES ────────────────────────────────────────────────────────

library(ggplot2)
library(dplyr)
library(gridExtra)
library(MASS)       # for polr() — ordinal logistic regression
library(pROC)

# Create images folder if it doesn't exist
dir.create("images", showWarnings = FALSE)

# ── 2. LOAD & CLEAN ──────────────────────────────────────────────────────────

ckd_raw <- read.csv(
  "C:/Users/Adharsh/OneDrive/Desktop/Desktop/Research/ckd-dataset-v2.csv",
  stringsAsFactors = FALSE,
  na.strings = c("", "?", "\t?")
)

# Drop first two junk rows (metadata + NA row)
ckd <- ckd_raw[3:nrow(ckd_raw), ]
rownames(ckd) <- NULL

cat("=== DATASET DIMENSIONS ===\n")
cat("Rows:", nrow(ckd), "| Columns:", ncol(ckd), "\n")
cat("Complete cases:", sum(complete.cases(ckd)), "\n\n")

# ── 3. DEFINE ORDERED FACTOR LEVELS ──────────────────────────────────────────

# Target
ckd$class <- factor(ckd$class, levels = c("notckd", "ckd"))

# CKD stage (ordinal: s1 = mildest, s5 = kidney failure)
ckd$stage <- factor(ckd$stage,
                    levels = c("s1", "s2", "s3", "s4", "s5"), ordered = TRUE)

# Hemoglobin (g/dL)
ckd$hemo <- factor(ckd$hemo,
                   levels = c("< 6.1", "6.1 - 7.4", "7.4 - 8.7", "8.7 - 10",
                              "10 - 11.3", "11.3 - 12.6", "12.6 - 13.9",
                              "13.9 - 15.2", "15.2 - 16.5", "≥ 16.5"),
                   ordered = TRUE)

# Serum Creatinine (mg/dL)
ckd$sc <- factor(ckd$sc,
                 levels = c("< 3.65", "3.65 - 7.29", "7.29 - 10.94",
                            "10.94 - 14.58", "≥ 14.58"),
                 ordered = TRUE)

# Blood Urea (mg/dL)
ckd$bu <- factor(ckd$bu,
                 levels = c("< 48.1", "48.1 - 86.2", "86.2 - 124.3",
                            "124.3 - 162.4", "≥ 162.4"),
                 ordered = TRUE)

# Albumin (ordinal 0-5 scale)
ckd$al <- factor(ckd$al,
                 levels = c("< 0", "1 - 1", "2 - 2", "3 - 3", "4 - 4", "≥ 4"),
                 ordered = TRUE)

# Blood Glucose Random (mg/dL)
ckd$bgr <- factor(ckd$bgr,
                  levels = c("< 112", "112 - 154", "154 - 196", "196 - 238",
                             "238 - 280", "≥ 280"),
                  ordered = TRUE)

# Age bins
ckd$age <- factor(ckd$age,
                  levels = c("< 12", "12 - 20", "20 - 27", "27 - 35", "35 - 43",
                             "43 - 51", "51 - 59", "59 - 66", "66 - 74", "≥ 74"),
                  ordered = TRUE)

# GFR (mL/min/1.73m²)
ckd$grf <- factor(ckd$grf,
                  levels = c("< 15", "15 - 29.1", "29.1 - 55.3", "55.3 - 89.8",
                             "89.8 - 127.3", "127.281 - 152.446",
                             "152.446 - 189.809", "189.809 - 227.944", "≥ 227.944"),
                  ordered = TRUE)

# Binary clinical flags
binary_cols <- c("rbc", "pc", "pcc", "ba", "htn", "dm", "cad",
                 "appet", "pe", "ane", "bp..Diastolic.", "bp.limit", "su")
for (col in binary_cols) ckd[[col]] <- as.factor(ckd[[col]])

cat("=== TARGET DISTRIBUTION ===\n")
print(table(ckd$class))
cat("CKD prevalence:", round(mean(ckd$class == "ckd") * 100, 1), "%\n\n")

cat("=== CKD STAGE DISTRIBUTION ===\n")
print(table(ckd$stage))
cat("\n")

# ── 4. EXPLORATORY ANALYSIS ───────────────────────────────────────────────────

# 4a. Stage distribution bar chart
stage_df <- as.data.frame(table(ckd$stage))
colnames(stage_df) <- c("stage", "count")
stage_df$label <- c("S1\nNormal/High\nGFR", "S2\nMildly\nDecreased",
                    "S3\nModerately\nDecreased", "S4\nSeverely\nDecreased",
                    "S5\nKidney\nFailure")

p_stage <- ggplot(stage_df, aes(x = stage, y = count, fill = stage)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
  geom_text(aes(label = count), vjust = -0.5, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("#2ecc71", "#f1c40f", "#e67e22", "#e74c3c", "#8e44ad")) +
  scale_x_discrete(labels = stage_df$label) +
  labs(title    = "Patient Distribution Across CKD Severity Stages",
       subtitle = "S1 = mildest (normal GFR) → S5 = kidney failure",
       x = "CKD Stage", y = "Number of Patients") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey50"))

ggsave("images/01_stage_distribution.png", p_stage, width = 10, height = 6, dpi = 150)
cat("Saved: images/01_stage_distribution.png\n")

# 4b. Hemoglobin vs CKD class
hemo_df <- ckd %>% group_by(hemo, class) %>% summarise(count = n(), .groups = "drop")

p_hemo <- ggplot(hemo_df, aes(x = hemo, y = count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("notckd" = "#2ecc71", "ckd" = "#e74c3c"),
                    labels = c("Non-CKD", "CKD")) +
  labs(title    = "Hemoglobin Level Distribution by CKD Status",
       subtitle = "Lower hemoglobin is a hallmark of renal anaemia in CKD",
       x = "Hemoglobin (g/dL)", y = "Count", fill = "Status") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 35, hjust = 1),
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey50"))

ggsave("images/02_hemoglobin_vs_class.png", p_hemo, width = 11, height = 6, dpi = 150)
cat("Saved: images/02_hemoglobin_vs_class.png\n")

# 4c. Serum Creatinine vs CKD class
sc_df <- ckd %>% group_by(sc, class) %>% summarise(count = n(), .groups = "drop")

p_sc <- ggplot(sc_df, aes(x = sc, y = count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("notckd" = "#2ecc71", "ckd" = "#e74c3c"),
                    labels = c("Non-CKD", "CKD")) +
  labs(title    = "Serum Creatinine Distribution by CKD Status",
       subtitle = "Rising creatinine reflects declining glomerular filtration rate",
       x = "Serum Creatinine (mg/dL)", y = "Count", fill = "Status") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 25, hjust = 1),
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey50"))

ggsave("images/03_creatinine_vs_class.png", p_sc, width = 10, height = 6, dpi = 150)
cat("Saved: images/03_creatinine_vs_class.png\n")

# 4d. Hemoglobin progression across CKD stages (among CKD patients only)
hemo_stage_df <- ckd %>%
  filter(class == "ckd") %>%
  group_by(stage, hemo) %>%
  summarise(count = n(), .groups = "drop")

p_hemo_stage <- ggplot(hemo_stage_df, aes(x = stage, y = hemo, size = count, colour = stage)) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c("s1" = "#2ecc71", "s2" = "#f1c40f",
                                 "s3" = "#e67e22", "s4" = "#e74c3c", "s5" = "#8e44ad")) +
  scale_size_continuous(range = c(3, 12)) +
  labs(title    = "Hemoglobin Levels Across CKD Severity Stages",
       subtitle = "Among CKD patients — bubble size = number of patients",
       x = "CKD Stage", y = "Hemoglobin (g/dL)",
       colour = "Stage", size = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey50"))

ggsave("images/04_hemo_across_stages.png", p_hemo_stage, width = 11, height = 7, dpi = 150)
cat("Saved: images/04_hemo_across_stages.png\n")

# 4e. GFR distribution across stages
grf_stage_df <- ckd %>% group_by(stage, grf) %>% summarise(count = n(), .groups = "drop")

p_grf <- ggplot(grf_stage_df, aes(x = stage, y = grf, size = count, colour = stage)) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c("s1" = "#2ecc71", "s2" = "#f1c40f",
                                 "s3" = "#e67e22", "s4" = "#e74c3c", "s5" = "#8e44ad")) +
  scale_size_continuous(range = c(3, 12)) +
  labs(title    = "GFR Distribution Across CKD Stages",
       subtitle = "GFR falls as stage worsens — confirming dataset consistency",
       x = "CKD Stage", y = "GFR (mL/min/1.73m²)",
       colour = "Stage", size = "Count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y  = element_text(size = 9),
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey50"))

ggsave("images/05_grf_across_stages.png", p_grf, width = 11, height = 7, dpi = 150)
cat("Saved: images/05_grf_across_stages.png\n")

# ── 5. CHI-SQUARE ASSOCIATION TESTS ──────────────────────────────────────────

cat("\n=== CHI-SQUARE TESTS: Association with CKD Class ===\n")
cat(sprintf("%-12s | %10s | %8s | %s\n", "Variable", "Chi-sq", "p-value", "Sig"))
cat(strrep("-", 50), "\n")

test_vars <- c("hemo", "sc", "bu", "al", "bgr", "htn", "dm", "ane", "pe", "grf")

for (v in test_vars) {
  tbl <- table(ckd[[v]], ckd$class)
  result <- tryCatch({
    cs <- suppressWarnings(chisq.test(tbl))
    if (is.nan(cs$statistic) || any(cs$expected < 5)) {
      ft <- fisher.test(tbl, simulate.p.value = TRUE, B = 10000)
      list(statistic = NA, p.value = ft$p.value, method = "Fisher")
    } else {
      list(statistic = cs$statistic, p.value = cs$p.value, method = "Chi-sq")
    }
  }, error = function(e) list(statistic = NA, p.value = NA, method = "Error"))
  sig      <- ifelse(is.na(result$p.value), "NA",
                     ifelse(result$p.value < 0.001, "***",
                            ifelse(result$p.value < 0.01,  "**",
                                   ifelse(result$p.value < 0.05,  "*", ""))))
  stat_str <- ifelse(is.na(result$statistic), "Fisher",
                     sprintf("%.2f", result$statistic))
  cat(sprintf("%-12s | %10s | %8.4f | %s\n", v, stat_str,
              ifelse(is.na(result$p.value), 0, result$p.value), sig))
}

# ── 6. ORDINAL LOGISTIC REGRESSION: Predicting CKD Stage ─────────────────────

cat("\n=== ORDINAL LOGISTIC REGRESSION: Predicting CKD Stage ===\n")
cat("(Among CKD patients — modelling severity progression s1 → s5)\n\n")

ckd_only <- ckd %>% filter(class == "ckd")

ord_model <- polr(
  stage ~ hemo + sc + bu + al,
  data   = ckd_only,
  Hess   = TRUE,
  method = "logistic"
)

print(summary(ord_model))

ctable <- coef(summary(ord_model))
p_vals <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
cat("\n--- Coefficients with p-values ---\n")
print(cbind(ctable, `p value` = round(p_vals, 4)))

# ── 7. BINARY LOGISTIC REGRESSION: CKD vs Non-CKD ────────────────────────────

cat("\n=== BINARY LOGISTIC REGRESSION: CKD vs Non-CKD ===\n")

# Convert ordered factors to numeric rank
ckd$hemo_n <- as.numeric(ckd$hemo)
ckd$sc_n   <- as.numeric(ckd$sc)
ckd$bu_n   <- as.numeric(ckd$bu)
ckd$al_n   <- as.numeric(ckd$al)
ckd$grf_n  <- as.numeric(ckd$grf)

# Use complete cases only to avoid ROC length mismatch
ckd_cc <- ckd[complete.cases(ckd[, c("class", "hemo_n", "bu_n", "al_n")]), ]
cat(sprintf("Complete cases for logistic model: %d\n", nrow(ckd_cc)))

bin_model <- glm(
  class ~ hemo_n + bu_n + al_n,
  data   = ckd_cc,
  family = binomial
)

print(summary(bin_model))

cat("\n--- Odds Ratios (with 95% CI) ---\n")
or_table <- exp(cbind(OR = coef(bin_model), confint(bin_model)))
print(round(or_table, 3))

# ── 8. ROC CURVE & AUC ────────────────────────────────────────────────────────

pred_probs <- predict(bin_model, type = "response")
roc_obj    <- roc(ckd_cc$class, pred_probs, levels = c("notckd", "ckd"), direction = "<")
auc_val    <- auc(roc_obj)

cat("\n=== AUC ===\n")
cat("AUC:", round(auc_val, 4), "\n")

png("images/06_roc_curve.png", width = 800, height = 700, res = 120)
plot(roc_obj,
     col = "#e74c3c", lwd = 2.5,
     main = paste0("ROC Curve — CKD Logistic Regression\nAUC = ", round(auc_val, 3)),
     xlab = "1 - Specificity (False Positive Rate)",
     ylab = "Sensitivity (True Positive Rate)",
     cex.main = 1.2, cex.lab = 1.1)
abline(a = 0, b = 1, col = "grey60", lty = 2, lwd = 1.5)
legend("bottomright",
       legend = c(paste0("Logistic Regression (AUC = ", round(auc_val, 3), ")"), "Random Chance"),
       col = c("#e74c3c", "grey60"), lty = c(1, 2), lwd = c(2.5, 1.5), cex = 0.9)
dev.off()
cat("Saved: images/06_roc_curve.png\n")

# Confusion matrix at 0.50
pred_class  <- ifelse(pred_probs >= 0.5, "ckd", "notckd")
cm          <- table(Predicted = pred_class, Actual = ckd_cc$class)
cat("\n=== CONFUSION MATRIX (threshold = 0.50) ===\n")
print(cm)

tp  <- cm["ckd", "ckd"];    tn <- cm["notckd", "notckd"]
fp  <- cm["ckd", "notckd"]; fn <- cm["notckd", "ckd"]
sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
accuracy    <- (tp + tn) / nrow(ckd_cc)

cat(sprintf("Accuracy    : %.1f%%\n", accuracy    * 100))
cat(sprintf("Sensitivity : %.1f%%\n", sensitivity * 100))
cat(sprintf("Specificity : %.1f%%\n", specificity * 100))

# ── 9. THRESHOLD ANALYSIS ─────────────────────────────────────────────────────

thresholds <- seq(0.1, 0.9, by = 0.05)
thresh_df  <- data.frame(threshold = thresholds, sensitivity = NA,
                         specificity = NA, accuracy = NA)

for (i in seq_along(thresholds)) {
  t    <- thresholds[i]
  prd  <- ifelse(pred_probs >= t, "ckd", "notckd")
  tp_  <- sum(prd == "ckd"    & ckd_cc$class == "ckd")
  tn_  <- sum(prd == "notckd" & ckd_cc$class == "notckd")
  fp_  <- sum(prd == "ckd"    & ckd_cc$class == "notckd")
  fn_  <- sum(prd == "notckd" & ckd_cc$class == "ckd")
  thresh_df$sensitivity[i] <- tp_ / (tp_ + fn_)
  thresh_df$specificity[i] <- tn_ / (tn_ + fp_)
  thresh_df$accuracy[i]    <- (tp_ + tn_) / nrow(ckd_cc)
}

thresh_long <- reshape(thresh_df,
                       varying = c("sensitivity", "specificity", "accuracy"),
                       v.names = "value", timevar = "metric",
                       times   = c("Sensitivity", "Specificity", "Accuracy"),
                       direction = "long")

p_thresh <- ggplot(thresh_long, aes(x = threshold, y = value, colour = metric)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = 0.35, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.38, y = 0.15, label = "Clinical\nscreening\nthreshold",
           size = 3.2, colour = "grey40") +
  scale_colour_manual(values = c("Sensitivity" = "#e74c3c",
                                 "Specificity" = "#2ecc71",
                                 "Accuracy"    = "#3498db")) +
  labs(title    = "Sensitivity vs. Specificity vs. Accuracy Across Thresholds",
       subtitle = "CKD Binary Logistic Regression",
       x = "Classification Threshold", y = "Metric Value", colour = "Metric") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey50"),
        legend.position = "bottom")

ggsave("images/07_threshold_tradeoff.png", p_thresh, width = 10, height = 6, dpi = 150)
cat("Saved: images/07_threshold_tradeoff.png\n")

# ── 10. SUMMARY ───────────────────────────────────────────────────────────────

cat("\n")
cat("=============================================================\n")
cat("                    ANALYSIS SUMMARY\n")
cat("=============================================================\n")
cat(sprintf("Total patients                   : %d\n", nrow(ckd)))
cat(sprintf("CKD patients                     : %d (%.1f%%)\n",
            sum(ckd$class == "ckd"), mean(ckd$class == "ckd") * 100))
cat(sprintf("AUC                              : %.3f\n", auc_val))
cat(sprintf("Accuracy    (threshold 0.50)     : %.1f%%\n", accuracy    * 100))
cat(sprintf("Sensitivity (threshold 0.50)     : %.1f%%\n", sensitivity * 100))
cat(sprintf("Specificity (threshold 0.50)     : %.1f%%\n", specificity * 100))
cat("=============================================================\n")
cat("\nKey clinical insight:\n")
cat("Hemoglobin declines progressively from Stage 1 to Stage 5,\n")
cat("reflecting worsening renal anaemia as erythropoietin production\n")
cat("falls. GFR and creatinine confirm the filtration collapse.\n")
cat("The staging variable in this dataset enables severity modelling\n")
cat("— a more clinically actionable question than CKD vs non-CKD.\n")
cat("=============================================================\n")