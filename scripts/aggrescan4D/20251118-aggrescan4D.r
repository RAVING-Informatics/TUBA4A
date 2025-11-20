#############################################
# TUBA4A AggreScan4D analysis
# - Reads residue-level AggreScan4D scores for WT and mutants
# - Computes average score per variant
# - Merges with variant labels (disease category)
# - Generates a lollipop plot vs WT aggregation
# - Exports a table with average scores (and a z-score column)
#############################################

## 1. Load packages ----
library(ggplot2)
library(dplyr)

## 2. Define input and output paths ----
# Adjust these to match your repository structure
input_dir   <- "data/aggrescan4d"
results_dir <- "results"
fig_dir     <- file.path(results_dir, "figures")
tab_dir     <- file.path(results_dir, "tables")

# Create output directories if they don't exist
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# Input files (tab-delimited text)
agg_file   <- file.path(input_dir, "TUBA4A-aggregation4D_scores.txt")
label_file <- file.path(input_dir, "variant-labels.txt")

# Output files
plot_file <- file.path(fig_dir, "Fig3_panel_G_aggrescan4D.png")
table_file <- file.path(tab_dir, "TUBA4A_aggrescan4D_avg_scores_z.txt")

## 3. Read data ----

# Residue-level AggreScan4D scores
# - First column = row names (residue index or identifier)
# - Remaining columns = variants
all_residue_scores <- read.delim(
  agg_file,
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Variant labels (variant + disease category)
variant_labels <- read.delim(
  label_file,
  header           = TRUE,
  stringsAsFactors = FALSE
)

## 4. Compute average score per variant ----

avg_vec <- colMeans(all_residue_scores, na.rm = TRUE)

avg_scores <- data.frame(
  variant   = names(avg_vec),
  Avg_score = as.numeric(avg_vec),
  stringsAsFactors = FALSE
)

## 5. Merge with labels ----

avg_scores <- left_join(avg_scores, variant_labels, by = "variant")

# Sanity check for unmatched variants
missing_labels <- dplyr::filter(avg_scores, is.na(label))
if (nrow(missing_labels) > 0) {
  warning("The following variants have no disease group label:")
  print(missing_labels)
}

## 6. Set ordering and factor levels ----

# Order of variants along x-axis (taken from label file)
variant_order <- variant_labels$variant

avg_scores$variant <- factor(avg_scores$variant, levels = variant_order)

avg_scores$label <- factor(
  avg_scores$label,
  levels = c("Wild-type",
             "Myopathy",
             "Multisystem proteinopathy",
             "ALS",
             "Ataxia",
             "Benign")
)

## 7. WT reference score and Y limits ----

wt_score <- avg_scores$Avg_score[avg_scores$variant == "WT"]

# Y-axis limits: slightly below min and slightly above WT
y_min <- min(avg_scores$Avg_score, na.rm = TRUE) - 0.01
y_max <- wt_score + 0.05

## 8. Define colour and shape mappings ----

avg_scores$Color <- ifelse(
  avg_scores$Avg_score > wt_score,
  "More aggregation than WT",
  "Less aggregation than WT"
)

variant_shapes <- c(
  "Wild-type"                  = 1,  # Empty circle
  "Myopathy"                   = 17, # Triangle
  "ALS"                        = 15, # Square
  "Ataxia"                     = 18, # Diamond
  "Multisystem proteinopathy"  = 10, # Crossed circle
  "Benign"                     = 4   # Cross
)

color_map <- c(
  "More aggregation than WT" = "red",
  "Less aggregation than WT" = "blue"
)

## 9. Create lollipop plot ----

plot_G <- ggplot(avg_scores, aes(x = variant, y = Avg_score)) +
  # Lollipop stems: from WT line to each point
  geom_segment(
    aes(xend = variant, y = wt_score, yend = Avg_score),
    colour = "black",
    linewidth = 0.65
  ) +
  # Points: coloured by aggregation vs WT, shaped by label (phenotype group)
  geom_point(aes(color = Color, shape = label), size = 3) +
  # Horizontal line at WT score
  geom_hline(
    yintercept = wt_score,
    linetype   = "dashed",
    color      = "green",
    linewidth  = 1
  ) +
  # Shaded regions for more vs less aggregation than WT
  annotate(
    "rect",
    xmin = 0.5,
    xmax = length(variant_order) + 0.5,
    ymin = wt_score,
    ymax = y_max,
    fill = "red",
    alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = 0.5,
    xmax = length(variant_order) + 0.5,
    ymin = y_min,
    ymax = wt_score,
    fill = "blue",
    alpha = 0.08
  ) +
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = variant_shapes) +
  scale_x_discrete(limits = variant_order) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(y_min, y_max, length.out = 5),
    labels = function(x) sprintf("%.2f", x)
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 12),
    legend.title = element_blank()
  ) +
  ylab("Average Aggrescan 4D score (normalised)") +
  ggtitle("Predicted aggregation propensity scores of WT and TUBA4A mutants")

# Show in R session
print(plot_G)

## 10. Save plot ----

ggsave(
  filename = plot_file,
  plot     = plot_G,
  width    = 12,
  height   = 8,
  dpi      = 300,
  units    = "in",
  bg       = "white"
)
cat("✓ Saved AggreScan4D plot to:", plot_file, "\n")