###############################################################
# TUBA4A – Genotype–Phenotype Correlation Analysis
# ------------------------------------------------------------
# This script computes correlations between TUBA4A variant
# domain location and clinical features.
#
# It:
#   1. Reads the curated clinical dataset text file
#   2. Coerces columns to the correct numeric/binary formats
#   3. Derives domain dummy variables and a limb Severity score
#   4. Computes Pearson correlations
#   5. Prints domain-wise summaries to the console
#   6. Produces a domain × clinical correlation heatmap
###############################################################

library(tidyverse)
library(ggplot2)

## ============================================================
## 1. IMPORT DATA
## ============================================================

input_file <- "data/genotype_phenotype/TUBA4A-clinical-correlation.txt"

data <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("", "NA", "n/a", "N/A")
)

## ------------------------------------------------------------
## IMPORTANT: Data Type Coercion
##
## MRC-based limb scores (DUL, PUL, DLL, PLL):
##   - Numeric, 1 = severe weakness, 5 = normal strength
##   - Higher = milder disease
##
## Binary clinical features (Axial, Facial, Ptosis, Respiratory):
##   - 1 = feature present, 0 = absent
##
## CK scale:
##   - 0 = normal
##   - 1 = mildly elevated
##   - 2 = markedly elevated
##
## Domain variable:
##   - Domain = 1 (GTPase), 2 (C-terminal), 3 (GTP-binding)
##   - In the analysis we convert this into domain-specific
##     dummy variables: Domain_GTPase, Domain_Cterminal,
##     Domain_GTPbinding (1 = variant in that domain, 0 = otherwise).
##
## Directionality matters:
##   - Limb MRC scores DECREASE with severity (1 = worst, 5 = best)
##   - Binary features and CK INCREASE with severity
## ------------------------------------------------------------

mrc_cols    <- c("DUL", "PUL", "DLL", "PLL")
binary_cols <- c("Axial", "Facial", "Ptosis", "Respiratory")
ck_col      <- "CK"

# Coerce all expected variables to numeric
data[mrc_cols]    <- lapply(data[mrc_cols],    function(x) as.numeric(as.character(x)))
data[binary_cols] <- lapply(data[binary_cols], function(x) as.numeric(as.character(x)))
data[[ck_col]]    <- as.numeric(as.character(data[[ck_col]]))

## ============================================================
## 2. DERIVE DOMAIN DUMMY VARIABLES
## ============================================================

data <- data %>%
  mutate(
    Domain_GTPase     = as.numeric(Domain == 1),
    Domain_Cterminal  = as.numeric(Domain == 2),
    Domain_GTPbinding = as.numeric(Domain == 3)
  )

## ============================================================
## 3. DEFINE SEVERITY SCORE
## ------------------------------------------------------------
## Severity is the average of the 4 limb MRC scores:
##
##   Severity = mean(DUL, PUL, DLL, PLL)
##
## Since all MRC scores range between 1–5:
##   - Severity also ranges between ~1–5
##   - HIGHER values = milder overall limb weakness
##   - LOWER values = more severe limb weakness
## ============================================================

data <- data %>%
  mutate(Severity = rowMeans(select(., all_of(mrc_cols)), na.rm = TRUE))

## ============================================================
## 4. BUILD CORRELATION DATASET
## ============================================================

cor_data <- data %>%
  select(
    # Limb MRC scores (continuous, reversed scale)
    all_of(mrc_cols),
    # Binary clinical features (0/1)
    all_of(binary_cols),
    # CK scale (0/1/2)
    CK,
    # Domain dummy variables
    Domain_GTPbinding,
    Domain_GTPase,
    Domain_Cterminal,
    # Derived overall severity score
    Severity
  ) %>%
  mutate(across(everything(), as.numeric))

## ============================================================
## 5. CORRELATION MATRIX
## ============================================================

cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

## For clarity, we focus on the biologically meaningful
## domain × clinical correlations (no domain–domain cells).

domain_vars <- c("Domain_GTPase", "Domain_Cterminal", "Domain_GTPbinding")
clinical_vars <- c(
  mrc_cols,
  "Severity",
  binary_cols,
  "CK"
)

domain_corr <- cor_matrix[domain_vars, clinical_vars, drop = FALSE]

cat("=== DOMAIN-CLINICAL CORRELATIONS (r, rounded to 2 decimals) ===\n")
print(round(domain_corr, 2))

## ============================================================
## 6. DOMAIN-WISE INTERPRETATION (CONSOLE + TEXT FILE)
## ============================================================

cat("\n=== DOMAIN-WISE INTERPRETATION ===\n\n")

summary_dir  <- "results/tables"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_file <- file.path(summary_dir, "geno_pheno_correlations_summary.txt")

# Capture printed output into a character vector
summary_lines <- c()

interpret_domain_to_text <- function(domain_name) {
  r_MRC <- cor_matrix[domain_name, mrc_cols, drop = FALSE]
  r_bin <- cor_matrix[domain_name, binary_cols, drop = FALSE]
  r_CK  <- cor_matrix[domain_name, "CK"]
  r_sev <- cor_matrix[domain_name, "Severity"]

  block <- c(
    "-------------------------------------------------",
    paste("Domain:", domain_name),
    "-------------------------------------------------",
    "",
    "LIMB WEAKNESS (MRC scores; higher = milder):",
    capture.output(print(round(r_MRC, 3))),
    "",
    "BINARY FEATURES (1 = present):",
    capture.output(print(round(r_bin, 3))),
    "",
    paste("CK correlation (0/1/2 scale): r =", round(r_CK, 3)),
    "",
    paste("Overall Severity (mean MRC score; higher = milder): r =", round(r_sev, 3)),
    ""
  )

  # print to console
  cat(paste(block, collapse = "\n"), "\n")

  return(block)
}

# Generate and collect the summaries
summary_lines <- c(
  "TUBA4A Clinical Correlation – Domain-wise Summary",
  "=================================================",
  ""
)

summary_lines <- c(
  summary_lines,
  interpret_domain_to_text("Domain_GTPase"),
  interpret_domain_to_text("Domain_Cterminal"),
  interpret_domain_to_text("Domain_GTPbinding")
)

# Write to file
writeLines(summary_lines, summary_file)
cat("✓ Saved domain-wise summary to:", summary_file, "\n")

## ============================================================
## 7. DOMAIN-CLINICAL FEATURE CORRELATION HEATMAP
## ============================================================

# Ensure output directory exists
results_dir <- "results"
fig_dir     <- file.path(results_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Reshape to long format for ggplot
domain_long <- as.data.frame(as.table(domain_corr))
colnames(domain_long) <- c("Domain", "Clinical_feature", "r")

# Nice readable domain names
domain_long$Domain <- factor(
  domain_long$Domain,
  levels = c("Domain_GTPase", "Domain_Cterminal", "Domain_GTPbinding"),
  labels = c("GTPase domain", "C-terminal domain", "GTP-binding region")
)

p_domain <- ggplot(domain_long, aes(x = Clinical_feature, y = Domain, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(
    low  = "#2166AC",
    mid  = "white",
    high = "#B2182B",
    limits = c(-1, 1),
    name = "r"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y  = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 14)
  ) +
  ggtitle("Genotype–Phenotype Correlation")

fig_file <- file.path(fig_dir, "Supp-Fig.correlation_genotype_vs_phenotype.png")

ggsave(
  filename = fig_file,
  plot     = p_domain,
  width    = 10,
  height   = 4.2,
  dpi      = 300,
  units    = "in",
  bg       = "white"
)

cat("✓ Saved heatmap to:", fig_file, "\n")
