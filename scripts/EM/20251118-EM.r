#############################################
# TUBA4A EM analysis – Figure 4B
# - Reads EM microtubule diameter measurements
# - Produces combined longitudinal + cross plot
# - Runs basic statistics (patients vs controls, per orientation)
#############################################

library(tidyverse)
library(ggdist)
library(patchwork)
library(effsize)   # Cohen's d
library(broom)     # tidy t-test / var-test output

## ===== 1. PATHS =====

input_dir   <- "data/em"
results_dir <- "results"
fig_dir     <- file.path(results_dir, "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

file <- file.path(input_dir, "TUBA4A-EM.txt")

## ===== 2. READ & CLEAN DATA (unchanged logic) =====

# Read everything as raw lines so we can skip the double header
raw_lines <- readLines(file)

# First two lines are headers, the rest is numeric data
# Write the numeric part into a temporary connection
tmp <- textConnection(raw_lines[-c(1, 2)])

em <- read.delim(
  tmp,
  header     = FALSE,
  sep        = "\t",
  na.strings = c("", "-", "NA")
)
close(tmp)

# Should be 8: P1_long, P1_cross, P2_long, P2_cross, C1_long, C1_cross, C2_long, C2_cross
colnames(em) <- c("P1_long",  "P1_cross",
                  "P2_long",  "P2_cross",
                  "C1_long",  "C1_cross",
                  "C2_long",  "C2_cross")

## ===== 3. TIDY TO LONG FORMAT (unchanged) =====

em_long <- em %>%
  mutate(row = row_number()) %>%   # just to keep row identity, if needed
  pivot_longer(
    cols      = -row,
    names_to  = c("Sample", "Orientation"),
    names_sep = "_",
    values_to = "Diameter_nm"
  ) %>%
  drop_na(Diameter_nm) %>%
  mutate(
    Sample = factor(
      Sample,
      levels = c("P1", "P2", "C1", "C2"),
      labels = c("Pat.9", "Pat.8", "Control.1", "Control.2")
    ),
    Orientation = factor(
      Orientation,
      levels = c("long", "cross"),
      labels = c("Longitudinal", "Cross")
    )
  )

## ===== 4. PLOTS PER INDIVIDUAL (ggdist) – unchanged =====

# Helper to avoid repetition
plot_em <- function(df, title_text) {
  ggplot(df, aes(x = Sample, y = Diameter_nm)) +
    # distribution (half-violin) + median & intervals
    ggdist::stat_halfeye(
      adjust         = 0.7,
      width          = 0.6,
      .width         = c(0.5, 0.95),   # inner = IQR, outer ≈ 95%
      slab_alpha     = 0.7,
      point_interval = median_qi
    ) +
    # raw data as dots on the other side
    ggdist::stat_dots(
      side          = "left",
      justification = 1.1,
      binwidth      = 0.2
    ) +
    theme_classic() +
    labs(
      x     = "",
      y     = "Microtubule diameter (nm)",
      title = title_text
    )
}

# (A) Longitudinal
p_long <- em_long %>%
  filter(Orientation == "Longitudinal") %>%
  plot_em("Longitudinal")

# (B) Cross
p_cross <- em_long %>%
  filter(Orientation == "Cross") %>%
  plot_em("Cross")

## ===== 5. Combine plots (patchwork) – unchanged layout =====

combined_plot <- (p_long | p_cross)

# ---- Save high-resolution output (only path changed) ----
fig_file <- file.path(fig_dir, "Figure4B_EM_combined.png")

ggsave(
  filename = fig_file,
  plot     = combined_plot,
  width    = 14,          # adjust as needed
  height   = 6,           # adjust as needed
  dpi      = 600,         # high res for journals
  units    = "in",
  bg       = "white"
)

combined_plot

## 5. Statistics used in main/supp text ----
# Patients: Pat.9 + Pat.8
# Controls: Control.1 + Control.2
# Orientation: longitudinal vs cross

EM_for_stats <- em_long %>%
  mutate(
    sample      = as.character(Sample),
    orientation = tolower(as.character(Orientation))
  ) %>%
  transmute(
    sample,
    orientation,
    diameter = Diameter_nm
  )

df_long  <- EM_for_stats %>% filter(orientation == "longitudinal")
df_cross <- EM_for_stats %>% filter(orientation == "cross")

patients_long  <- df_long  %>% filter(grepl("Pat.9|Pat.8", sample))
controls_long  <- df_long  %>% filter(grepl("Control", sample))

patients_cross <- df_cross %>% filter(grepl("Pat.9|Pat.8", sample))
controls_cross <- df_cross %>% filter(grepl("Control", sample))

# Helper to compute t-test, Cohen's d, variance test
run_stats <- function(group1, group2, orientation_label) {
  t_res <- t.test(group1$diameter, group2$diameter, var.equal = FALSE)
  d_res <- cohen.d(group1$diameter, group2$diameter)
  f_res <- var.test(group1$diameter, group2$diameter)
  
  tibble(
    orientation      = orientation_label,
    mean_patients    = mean(group1$diameter, na.rm = TRUE),
    mean_controls    = mean(group2$diameter, na.rm = TRUE),
    t_statistic      = t_res$statistic,
    df               = t_res$parameter,
    p_value          = t_res$p.value,
    ci_lower         = t_res$conf.int[1],
    ci_upper         = t_res$conf.int[2],
    cohen_d          = as.numeric(d_res$estimate),
    cohen_d_magnitude = d_res$magnitude,
    var_ratio        = f_res$estimate,
    var_test_pvalue  = f_res$p.value
  )
}

stats_long  <- run_stats(patients_long,  controls_long,  "longitudinal")
stats_cross <- run_stats(patients_cross, controls_cross, "cross")

EM_stats <- bind_rows(stats_long, stats_cross)

# Inspect in R session
EM_stats

# Save stats for manuscript traceability
write.table(
  EM_stats,
  file      = stats_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)