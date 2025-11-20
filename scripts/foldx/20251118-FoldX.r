#############################################
# TUBA4A FoldX analysis
# - MYOPATHY VARIANTS:
#     * read FoldX ΔΔG values
#     * generate bar plot with error bars and pLDDT-based colours
#     * compute z-scores (stability_z) within this set
#     * export table with z-scores
# - NON-MYOPATHY VARIANTS:
#     * read FoldX ΔΔG values
#     * compute z-scores (stability_z) within this set
#     * export table with z-scores (no plot)
#############################################

## 1. Load packages ----
library(ggplot2)

## 2. Define input and output paths ----
# Adjust these to match your repo layout
input_dir   <- "data/foldx"
results_dir <- "results"
fig_dir     <- file.path(results_dir, "figures")
tab_dir     <- file.path(results_dir, "tables")

# Create output directories if they don't exist
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# Input file names (tab-delimited)
file_myopathy    <- file.path(input_dir, "TUBA4A_FoldX_myopathy.txt")
file_nonmyopathy <- file.path(input_dir, "TUBA4A_FoldX_nonmyopathy.txt")

# Output files
plot_file_myo    <- file.path(fig_dir, "Fig3_panel_F_foldx_myopathy.png")
table_myo_file   <- file.path(tab_dir, "TUBA4A_FoldX_zscore_myopathy.txt")
table_non_file   <- file.path(tab_dir, "TUBA4A_FoldX_zscore_nonmyopathy.txt")

#############################################
## 3. MYOPATHY VARIANTS: plot + z-scores ----
#############################################

# 3.1 Import data
foldX_myo <- read.delim(
  file_myopathy,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 3.2 Basic coercions
foldX_myo$diff  <- as.numeric(foldX_myo$diff)
foldX_myo$SD    <- as.numeric(foldX_myo$SD)
foldX_myo$pLDDT <- as.numeric(foldX_myo$pLDDT)

# Preserve order of AA on x-axis (unique to avoid duplicated-level errors)
foldX_myo$AA <- factor(foldX_myo$AA, levels = unique(foldX_myo$AA))

# 3.3 Define ΔΔG threshold
threshold <- 1.6  # kcal/mol, conservative threshold used in the manuscript

# 3.4 Plot: ΔΔG for myopathy variants
plot_F_myo <- ggplot(
  foldX_myo,
  aes(
    x    = AA,
    y    = diff,
    fill = factor(ifelse(pLDDT >= 90, "Very_confident", "Confident"))
  )
) +
  geom_col(show.legend = FALSE, colour = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept =  threshold,  linetype = "dashed", color = "red",
             linewidth = 0.8, alpha = 0.5) +
  geom_hline(yintercept = -threshold,  linetype = "dashed", color = "red",
             linewidth = 0.8, alpha = 0.5) +
  scale_x_discrete() +
  scale_y_continuous(breaks = seq(-6, 12, by = 2)) +
  scale_fill_manual(
    name   = "pLDDT",
    values = c("Confident" = "steelblue1", "Very_confident" = "dodgerblue4")
  ) +
  geom_errorbar(
    aes(ymin = diff - SD, ymax = diff + SD),
    width    = 0.4,
    colour   = "black",
    alpha    = 0.9,
    linewidth = 0.8
  ) +
  labs(
    y     = "Free energy change (kcal/mol)",
    title = "Predicted free energy changes for missense variants on protein stability"
  ) +
  theme_classic() +
  theme(
    axis.line        = element_line(linewidth = 1),
    axis.text.x      = element_text(angle = 60, hjust = 1),
    axis.text        = element_text(colour = "black", size = 14),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(face = "bold", size = 17),
    axis.ticks       = element_line(linewidth = 0.7),
    axis.ticks.length = unit(5, "pt"),
    plot.title       = element_text(size = 16, face = "bold", hjust = 0.5)
  )

print(plot_F_myo)

# 3.5 Save plot
ggsave(
  filename = plot_file_myo,
  plot     = plot_F_myo,
  width    = 12,
  height   = 8,
  dpi      = 300,
  units    = "in",
  bg       = "white"
)
cat("✓ Saved FoldX myopathy plot to:", plot_file_myo, "\n")

# 3.6 Compute z-scores for myopathy variants
# Z-scores are computed within this set only (no mixing with non-myopathy)
foldX_myo$stability_z <- as.numeric(scale(foldX_myo$diff))

# 3.7 Export myopathy table with z-scores
write.table(
  foldX_myo,
  file      = table_myo_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

#############################################
## 4. NON-MYOPATHY VARIANTS: z-scores only ----
#############################################

# 4.1 Import data
foldX_non <- read.delim(
  file_nonmyopathy,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 4.2 Basic coercions
foldX_non$diff  <- as.numeric(foldX_non$diff)
foldX_non$SD    <- as.numeric(foldX_non$SD)
foldX_non$pLDDT <- as.numeric(foldX_non$pLDDT)

# 4.3 Compute z-scores for non-myopathy variants
foldX_non$stability_z <- as.numeric(scale(foldX_non$diff))

# 4.4 Export non-myopathy table with z-scores
write.table(
  foldX_non,
  file      = table_non_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)