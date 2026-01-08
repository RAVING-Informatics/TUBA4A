###############################################################
# TUBA4A – Mechanistic Correlation Analysis
# ------------------------------------------------------------
# This script relates mechanistic readouts (FoldX ΔΔG z-scores,
# Aggrescan4D aggregation, and EM morphology) to TUBA4A domain
# location. It:
#   1. Reads the curated mechanistic dataset
#   2. Coerces variables to numeric and derives morphology severity
#   3. Defines two subsets:
#        Set A – variants with ΔΔG and aggregation scores only and no morphology data
#        Set B – variants with morphology readouts along with Set A
#   4. Computes Spearman correlations for both sets
#   5. Computes partial correlations (mechanistic vs morphology
#      controlling for domain)
#   6. Prints summaries to console and to a text file
#   7. Produces a combined Set A + Set B correlation heatmap
###############################################################

## ===== PACKAGES =====
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(corrplot)

## ============================================================
## 1. PATHS AND INPUT
## ============================================================

input_file  <- "data/mechanistic_correlation/TUBA4A-mechanistic-correlation_updated.txt"
results_dir <- "results"
fig_dir     <- file.path(results_dir, "figures")
tab_dir     <- file.path(results_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

stats_file  <- file.path(tab_dir, "mechanistic_correlations_summary.txt")
fig_file    <- file.path(fig_dir, "Supp-Fig.mechanistic_correlations.png")

## ============================================================
## 2. READ & RENAME
## ============================================================

raw <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("", "NA", "n/a", "N/A", "na")
)

# Expected columns (names in the .txt file are kept explicit here)
data <- raw %>%
  rename(
    Variant    = Variant,
    Domain     = Domain,   # 1 = GTPase, 2 = C-terminal, 3 = GTP-binding
    z_ddg      = Free.energy.change.scaled.z.score,
    agg_prop   = Aggregation.propensity.for.mutant.protein,
    pct_normal = pct_normal.microtubules,
    pct_mild   = pct_mild.microtubules,
    pct_abn    = pct_abnormal.microtubules
  )

# Coerce numeric variables
num_cols <- c("Domain", "z_ddg", "agg_prop",
              "pct_normal", "pct_mild", "pct_abn")

data[num_cols] <- lapply(data[num_cols], function(x) suppressWarnings(as.numeric(x)))

## ------------------------------------------------------------
## DOMAIN FACTOR & DUMMY VARIABLES
##
## Domain coding in the raw file:
##   1 = GTPase domain
##   2 = C-terminal domain
##   3 = GTP-binding region
##
## For convenience in analyses, we derive:
##   Domain_factor : factor("GTPase", "C_terminal", "GTP_binding")
##   D_GTPase, D_Cterminal, D_GTPbinding : 0/1 indicators
## ------------------------------------------------------------

data <- data %>%
  mutate(
    Domain_factor = factor(
      Domain,
      levels = c(1, 2, 3),
      labels = c("GTPase", "C_terminal", "GTP_binding")
    ),
    D_GTPase     = as.numeric(Domain_factor == "GTPase"),
    D_Cterminal  = as.numeric(Domain_factor == "C_terminal"),
    D_GTPbinding = as.numeric(Domain_factor == "GTP_binding")
  )

## ------------------------------------------------------------
## MORPHOLOGY SEVERITY SCORE
##
## pct_normal, pct_mild, pct_abn represent the percentage of
## microtubules classified as normal / mildly abnormal / abnormal.
##
## We build a composite morphology severity score on a 0–2 scale:
##   morph_severity =
##      (0 * %normal + 1 * %mild + 2 * %abnormal) / 100
##
## Interpretation:
##   0   = entirely normal morphology
##   1   = on average "mild"
##   2   = entirely abnormal morphology
## ------------------------------------------------------------

data <- data %>%
  mutate(
    morph_severity = (0 * pct_normal + 1 * pct_mild + 2 * pct_abn) / 100
  )

## ============================================================
## 3. DEFINE SUBSETS
## ============================================================

# Set A: variants with both ΔΔG and aggregation available
setA <- data %>%
  filter(!is.na(z_ddg) & !is.na(agg_prop))

# Set B: variants with morphology data available
setB <- data %>%
  filter(!is.na(morph_severity))

## ============================================================
## 4. HELPER: SAFE SPEARMAN CORRELATION
## ============================================================

spearman_pair <- function(df, x, y) {
  dd <- df %>% select(all_of(c(x, y))) %>% tidyr::drop_na()
  if (nrow(dd) < 3) return(NA_real_)  # too few points to be meaningful
  suppressWarnings(cor(dd[[1]], dd[[2]], method = "spearman"))
}

## ============================================================
## 5. SET A – DOMAIN, ΔΔG, AGGREGATION
## ============================================================

cat("\n=== SET A (ΔΔG & aggregation present) ===\n")

# Domain vs ΔΔG (treat domain as 1<2<3 for monotonic trend)
rho_A_dom_ddg <- spearman_pair(setA, "Domain", "z_ddg")
# Domain vs aggregation
rho_A_dom_agg <- spearman_pair(setA, "Domain", "agg_prop")

cat(sprintf("• Domain vs ΔΔG (z)             : ρ = %.3f\n", rho_A_dom_ddg))
cat(sprintf("• Domain vs aggregation         : ρ = %.3f\n", rho_A_dom_agg))

# Correlation matrix for Set A
corA_vars <- c("z_ddg", "agg_prop", "Domain")
corA <- setA %>% select(all_of(corA_vars))

# Remove zero-variance columns (just in case)
keepA <- sapply(corA, function(x) {
  xx <- x[!is.na(x)]
  length(xx) > 1 && var(xx) > 0
})
corA  <- corA[, keepA, drop = FALSE]

cor_mat_A <- cor(corA, method = "spearman", use = "pairwise.complete.obs")
cat("\nCorrelation matrix – Set A (Spearman ρ):\n")
print(round(cor_mat_A, 2))

## ============================================================
## 6. SET B – INCLUDING MORPHOLOGY
## ============================================================

cat("\n=== SET B (morphology present) ===\n")

rho_B_dom_ddg  <- spearman_pair(setB, "Domain", "z_ddg")
rho_B_dom_agg  <- spearman_pair(setB, "Domain", "agg_prop")
rho_B_dom_mph  <- spearman_pair(setB, "Domain", "morph_severity")
rho_B_ddg_mph  <- spearman_pair(setB, "z_ddg",  "morph_severity")
rho_B_agg_mph  <- spearman_pair(setB, "agg_prop","morph_severity")

cat(sprintf("• Domain vs ΔΔG (z)             : ρ = %.3f\n", rho_B_dom_ddg))
cat(sprintf("• Domain vs aggregation         : ρ = %.3f\n", rho_B_dom_agg))
cat(sprintf("• Domain vs morphology severity : ρ = %.3f\n", rho_B_dom_mph))
cat(sprintf("• ΔΔG (z) vs morphology severity: ρ = %.3f\n", rho_B_ddg_mph))
cat(sprintf("• Aggregation vs morphology sev.: ρ = %.3f\n", rho_B_agg_mph))

corB_vars <- c("z_ddg", "agg_prop", "morph_severity", "Domain")
corB <- setB %>% select(all_of(corB_vars))

keepB <- sapply(corB, function(x) {
  xx <- x[!is.na(x)]
  length(xx) > 1 && var(xx) > 0
})
corB  <- corB[, keepB, drop = FALSE]

cor_mat_B <- cor(corB, method = "spearman", use = "pairwise.complete.obs")
cat("\nCorrelation matrix – Set B (Spearman ρ):\n")
print(round(cor_mat_B, 2))

## ============================================================
## 7. PARTIAL CORRELATIONS (CONTROLLING FOR DOMAIN)
## ============================================================

cat("\n=== PARTIAL CORRELATIONS (| mechanistic vs morphology, controlling for domain) ===\n")

rho_pc1 <- NA_real_
rho_pc2 <- NA_real_

# ΔΔG vs morphology severity | domain
df_pc1 <- setB %>% select(z_ddg, morph_severity, Domain) %>% drop_na()
if (nrow(df_pc1) >= 3) {
  m1_z   <- lm(z_ddg ~ factor(Domain), data = df_pc1)
  m1_mph <- lm(morph_severity ~ factor(Domain), data = df_pc1)
  rho_pc1 <- suppressWarnings(cor(resid(m1_z), resid(m1_mph), method = "spearman"))
  cat(sprintf("• Partial ρ (ΔΔG vs morphology | domain)        : ρ = %.3f\n", rho_pc1))
} else {
  cat("• Partial ρ (ΔΔG vs morphology | domain)        : not enough data\n")
}

# Aggregation vs morphology severity | domain
df_pc2 <- setB %>% select(agg_prop, morph_severity, Domain) %>% drop_na()
if (nrow(df_pc2) >= 3) {
  m2_a   <- lm(agg_prop ~ factor(Domain), data = df_pc2)
  m2_mph <- lm(morph_severity ~ factor(Domain), data = df_pc2)
  rho_pc2 <- suppressWarnings(cor(resid(m2_a), resid(m2_mph), method = "spearman"))
  cat(sprintf("• Partial ρ (Aggregation vs morphology | domain): ρ = %.3f\n", rho_pc2))
} else {
  cat("• Partial ρ (Aggregation vs morphology | domain): not enough data\n")
}

## ============================================================
## 8. SAVE SUMMARY STATS TO TEXT FILE
## ============================================================

summary_lines <- c(
  "TUBA4A mechanistic correlations (Spearman ρ)",
  "=========================================",
  "",
  "SET A (ΔΔG & aggregation present)",
  sprintf("Domain vs ΔΔG (z)             : %.3f", rho_A_dom_ddg),
  sprintf("Domain vs aggregation         : %.3f", rho_A_dom_agg),
  "",
  "Correlation matrix – Set A:",
  capture.output(print(round(cor_mat_A, 3))),
  "",
  "SET B (morphology present)",
  sprintf("Domain vs ΔΔG (z)             : %.3f", rho_B_dom_ddg),
  sprintf("Domain vs aggregation         : %.3f", rho_B_dom_agg),
  sprintf("Domain vs morphology severity : %.3f", rho_B_dom_mph),
  sprintf("ΔΔG (z) vs morphology severity: %.3f", rho_B_ddg_mph),
  sprintf("Aggregation vs morphology sev.: %.3f", rho_B_agg_mph),
  "",
  "Correlation matrix – Set B:",
  capture.output(print(round(cor_mat_B, 3))),
  "",
  "PARTIAL CORRELATIONS (controlling for domain)",
  sprintf("ΔΔG vs morphology | domain        : %s",
          ifelse(is.na(rho_pc1), "NA", sprintf("%.3f", rho_pc1))),
  sprintf("Aggregation vs morphology | domain: %s",
          ifelse(is.na(rho_pc2), "NA", sprintf("%.3f", rho_pc2)))
)

writeLines(summary_lines, con = stats_file)
cat("\n✓ Saved correlation summary to:", stats_file, "\n")

## ============================================================
## 9. COMBINED HEATMAP (SET A + SET B)
## ============================================================
# We plot the two correlation matrices side-by-side using corrplot.

png(fig_file, width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))

corrplot(
  cor_mat_A,
  method = "color",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.9,
  addCoef.col = "black",
  number.cex = 0.7,
  col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
  title = "Set A: Domain, ΔΔG (z), Aggregation",
  mar = c(0, 0, 2, 0)
)

corrplot(
  cor_mat_B,
  method = "color",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.9,
  addCoef.col = "black",
  number.cex = 0.7,
  col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
  title = "Set B: Domain, ΔΔG (z), Aggregation, Morphology",
  mar = c(0, 0, 2, 0)
)

dev.off()
par(mfrow = c(1, 1))  # reset

cat("✓ Saved combined heatmap to:", fig_file, "\n")
cat("\n✓ Mechanistic correlation analysis complete.\n")