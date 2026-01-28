# =============================================================================
# PLOTTING SCRIPT FOR INTEGRATED POPULATION GENOMICS ANALYSIS
# =============================================================================
# This script creates publication-ready figures from the CSV outputs of the
# integrated population genomics analysis (ROH, Private Alleles, FST).
#
# Use this script when:
# - You've already run the main analysis and have the CSV files
# - You want to recreate figures with different parameters/colors
# - You want to customize figure layouts without re-running the analysis
#
# Author: Created for Michael Kantar
# Date: 2026-01-27
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# =============================================================================
# SECTION 1: PARAMETERS - CUSTOMIZE THESE
# =============================================================================

cat("=============================================================================\n")
cat("CREATING FIGURES FROM INTEGRATED ANALYSIS OUTPUTS\n")
cat("=============================================================================\n\n")

# --- Input CSV files (from the integrated analysis) ---
froh_file <- "froh_coefficients.csv"
private_allele_file <- "private_allele_summary.csv"
combined_results_file <- "fst_private_allele_analysis_complete.csv"

# --- Figure output settings ---
output_prefix <- "publication"  # Prefix for output files
figure_width <- 14              # inches
figure_height <- 10             # inches
figure_dpi <- 300               # DPI for standard version
figure_dpi_hires <- 600         # DPI for high-resolution version

# --- FST outlier threshold (must match what was used in analysis) ---
fst_threshold_995 <- 0.995
fst_threshold_99 <- 0.99

# --- Color schemes ---
# Panel A: Private alleles
private_allele_colors <- "viridis"  # Options: "viridis", "magma", "plasma", "inferno", "cividis"

# Panel B: Froh
froh_colors <- "plasma"  # Options: "viridis", "magma", "plasma", "inferno", "cividis"

# Panel C: Manhattan plot
manhattan_colors <- list(
  neutral = "gray60",
  private_only = "#2CA02C",      # Green
  fst_only = "#FF7F0E",          # Orange
  both = "#D62728"               # Red
)

# --- Text sizes ---
base_text_size <- 10
title_text_size <- 11
axis_text_size <- 9
axis_text_angle <- 45  # For x-axis labels

# --- Additional plot options ---
show_individual_points_froh <- TRUE  # Show individual data points on Froh boxplot
show_counts_on_bars <- TRUE          # Show counts above private allele bars
outlier_line_size <- 0.5             # Thickness of FST threshold lines

# =============================================================================
# SECTION 2: LOAD DATA
# =============================================================================

cat("Loading data files...\n")

# Load Froh data
if (!file.exists(froh_file)) {
  stop(sprintf("Error: Cannot find %s\nMake sure you've run the integrated analysis first!", froh_file))
}
froh_data <- read.csv(froh_file, stringsAsFactors = FALSE)
cat(sprintf("  ✓ Froh data loaded: %d individuals\n", nrow(froh_data)))

# Load private allele data
if (!file.exists(private_allele_file)) {
  stop(sprintf("Error: Cannot find %s\nMake sure you've run the integrated analysis first!", private_allele_file))
}
private_data <- read.csv(private_allele_file, stringsAsFactors = FALSE)
cat(sprintf("  ✓ Private allele data loaded: %d populations\n", nrow(private_data)))

# Load combined FST and private allele results
if (!file.exists(combined_results_file)) {
  stop(sprintf("Error: Cannot find %s\nMake sure you've run the integrated analysis first!", combined_results_file))
}
combined_data <- read.csv(combined_results_file, stringsAsFactors = FALSE)
cat(sprintf("  ✓ Combined FST/Private allele data loaded: %d SNPs\n", nrow(combined_data)))

cat("\nData summary:\n")
cat(sprintf("  Populations: %d\n", length(unique(froh_data$population))))
cat(sprintf("  Individuals: %d\n", nrow(froh_data)))
cat(sprintf("  Total private alleles: %d\n", sum(private_data$private_alleles)))
cat(sprintf("  SNPs analyzed: %d\n", nrow(combined_data)))
cat(sprintf("  FST outliers (99.5%%): %d\n", sum(combined_data$Outlier_995, na.rm = TRUE)))
cat(sprintf("  Private alleles in FST data: %d\n\n", sum(combined_data$Is_Private_Any, na.rm = TRUE)))

# =============================================================================
# SECTION 3: CREATE PANEL A - PRIVATE ALLELES BY POPULATION
# =============================================================================

cat("Creating Panel A: Private alleles by population...\n")

panel_a <- ggplot(private_data, 
                  aes(x = reorder(population, -private_alleles), 
                      y = private_alleles, 
                      fill = population)) +
  geom_col(color = "black", size = 0.3) +
  labs(
    title = "A) Private Alleles by Population",
    x = "Population",
    y = "Number of Private Alleles"
  ) +
  scale_fill_viridis_d(option = private_allele_colors) +
  theme_classic(base_size = base_text_size) +
  theme(
    axis.text.x = element_text(angle = axis_text_angle, hjust = 1, size = axis_text_size),
    axis.text.y = element_text(size = axis_text_size),
    axis.title = element_text(size = base_text_size, face = "bold"),
    plot.title = element_text(size = title_text_size, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )

# Optionally add count labels on bars
if (show_counts_on_bars) {
  panel_a <- panel_a +
    geom_text(aes(label = private_alleles), 
              vjust = -0.5, size = 3)
}

cat("  ✓ Panel A created\n")

# =============================================================================
# SECTION 4: CREATE PANEL B - FROH BY POPULATION
# =============================================================================

cat("Creating Panel B: Froh by population...\n")

panel_b <- ggplot(froh_data, 
                  aes(x = reorder(population, -Froh, FUN = median), 
                      y = Froh, 
                      fill = population)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               fill = "white", color = "black") +
  labs(
    title = "B) Inbreeding Coefficient (Froh) by Population",
    x = "Population",
    y = expression(F[ROH])
  ) +
  scale_fill_viridis_d(option = froh_colors) +
  theme_classic(base_size = base_text_size) +
  theme(
    axis.text.x = element_text(angle = axis_text_angle, hjust = 1, size = axis_text_size),
    axis.text.y = element_text(size = axis_text_size),
    axis.title = element_text(size = base_text_size, face = "bold"),
    plot.title = element_text(size = title_text_size, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )

# Optionally add individual points
if (show_individual_points_froh) {
  panel_b <- panel_b +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8)
}

cat("  ✓ Panel B created\n")

# =============================================================================
# SECTION 5: CREATE PANEL C - MANHATTAN PLOT
# =============================================================================

cat("Creating Panel C: Manhattan plot of FST with private alleles...\n")

# Calculate cumulative positions for Manhattan plot
chr_lengths <- combined_data %>%
  group_by(CHROM) %>%
  summarise(chr_length = max(POS), .groups = "drop") %>%
  arrange(CHROM) %>%
  mutate(
    chr_start = cumsum(lag(chr_length, default = 0)),
    chr_mid = chr_start + chr_length / 2,
    chr_end = chr_start + chr_length
  )

# Add cumulative positions to data
manhattan_data <- combined_data %>%
  left_join(chr_lengths %>% select(CHROM, chr_start), by = "CHROM") %>%
  mutate(pos_cumul = chr_start + POS) %>%
  mutate(
    point_type = case_when(
      Outlier_995 & Is_Private_Any ~ "FST outlier + Private",
      Outlier_995 & !Is_Private_Any ~ "FST outlier only",
      !Outlier_995 & Is_Private_Any ~ "Private allele only",
      TRUE ~ "Neutral"
    ),
    point_type = factor(point_type, 
                        levels = c("Neutral", "Private allele only", 
                                   "FST outlier only", "FST outlier + Private"))
  )

# Calculate actual thresholds from the data
threshold_995 <- quantile(combined_data$FST, probs = fst_threshold_995, na.rm = TRUE)
threshold_99 <- quantile(combined_data$FST, probs = fst_threshold_99, na.rm = TRUE)
threshold_95 <- quantile(combined_data$FST, probs = 0.95, na.rm = TRUE)

cat(sprintf("  FST threshold (99.5%%): %.4f\n", threshold_995))
cat(sprintf("  FST threshold (99%%): %.4f\n", threshold_99))
cat(sprintf("  FST threshold (95%%): %.4f\n", threshold_95))

# Create Manhattan plot - matching uploaded figure style
panel_c <- ggplot(manhattan_data, aes(x = pos_cumul / 1e6, y = FST)) +
  # All SNPs as background (light blue/gray)
  geom_point(aes(color = "Background"), size = 0.8, alpha = 0.4) +
  # Private alleles as GREEN TRIANGLES (overlaid on top)
  geom_point(data = manhattan_data %>% filter(Is_Private_Any == TRUE),
             aes(color = "Private alleles"),
             shape = 17,  # Triangle
             size = 1.5, 
             alpha = 0.8) +
  # Threshold lines
  geom_hline(yintercept = threshold_995, linetype = "dashed", 
             color = "red", size = outlier_line_size,
             alpha = 0.8) +
  geom_hline(yintercept = threshold_99, linetype = "dashed", 
             color = "orange", size = outlier_line_size,
             alpha = 0.8) +
  geom_hline(yintercept = threshold_95, linetype = "dashed", 
             color = "yellow3", size = outlier_line_size,
             alpha = 0.8) +
  # Color scale - simple two-color scheme
  scale_color_manual(
    values = c(
      "Background" = "lightblue3",
      "Private alleles" = "darkgreen"
    ),
    name = NULL,
    labels = c(
      "Background" = "All SNPs",
      "Private alleles" = "Private alleles"
    )
  ) +
  # Chromosome labels
  scale_x_continuous(
    breaks = chr_lengths$chr_mid / 1e6,
    labels = chr_lengths$CHROM,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "C) FST Manhattan Plot with Private Alleles",
    subtitle = "Green triangles = Private alleles | Red/Orange/Yellow = 99.5%/99%/95% FST thresholds",
    x = "Chromosome",
    y = "FST"
  ) +
  theme_classic(base_size = base_text_size) +
  theme(
    axis.text.x = element_text(angle = axis_text_angle, hjust = 1, size = axis_text_size - 1),
    axis.text.y = element_text(size = axis_text_size),
    axis.title = element_text(size = base_text_size, face = "bold"),
    plot.title = element_text(size = title_text_size, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = axis_text_size, hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = axis_text_size),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.border = element_rect(fill = NA, color = "gray50"),
    panel.background = element_rect(fill = "white")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, shape = c(16, 17))))

cat("  ✓ Panel C created\n")

# =============================================================================
# SECTION 6: COMBINE PANELS AND SAVE FIGURE
# =============================================================================

cat("\nCombining panels into final figure...\n")

combined_figure <- (panel_a | panel_b) / panel_c +
  plot_layout(heights = c(1, 1.2))

# Save standard resolution
output_file <- paste0(output_prefix, "_integrated_figure.png")
ggsave(output_file, 
       combined_figure, 
       width = figure_width, 
       height = figure_height, 
       dpi = figure_dpi)
cat(sprintf("  ✓ Saved: %s (%d x %d inches, %d DPI)\n", 
            output_file, figure_width, figure_height, figure_dpi))

# Save high-resolution version
output_file_hires <- paste0(output_prefix, "_integrated_figure_hires.png")
ggsave(output_file_hires, 
       combined_figure, 
       width = figure_width, 
       height = figure_height, 
       dpi = figure_dpi_hires)
cat(sprintf("  ✓ Saved: %s (%d x %d inches, %d DPI)\n\n", 
            output_file_hires, figure_width, figure_height, figure_dpi_hires))

# =============================================================================
# SECTION 7: CREATE INDIVIDUAL PANEL FILES (OPTIONAL)
# =============================================================================

cat("Saving individual panels...\n")

# Panel A only
ggsave(paste0(output_prefix, "_panel_A_private_alleles.png"), 
       panel_a, width = 7, height = 5, dpi = figure_dpi)

# Panel B only
ggsave(paste0(output_prefix, "_panel_B_froh.png"), 
       panel_b, width = 7, height = 5, dpi = figure_dpi)

# Panel C only
ggsave(paste0(output_prefix, "_panel_C_manhattan.png"), 
       panel_c, width = 14, height = 6, dpi = figure_dpi)

cat("  ✓ Individual panels saved\n\n")

# =============================================================================
# SECTION 8: CREATE SUPPLEMENTARY FIGURES
# =============================================================================

cat("Creating supplementary figures...\n")

# --- Supplementary Figure 1: Froh distribution histogram ---
supp_fig1 <- ggplot(froh_data, aes(x = Froh, fill = population)) +
  geom_histogram(bins = 30, alpha = 0.7, color = "black", size = 0.2) +
  facet_wrap(~population, scales = "free_y", ncol = 3) +
  labs(
    title = "Distribution of Froh by Population",
    x = expression(F[ROH]),
    y = "Count"
  ) +
  scale_fill_viridis_d(option = froh_colors) +
  theme_classic(base_size = base_text_size) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = axis_text_size),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )

ggsave(paste0(output_prefix, "_supplementary_froh_distribution.png"), 
       supp_fig1, width = 12, height = 8, dpi = figure_dpi)
cat("  ✓ Supplementary Figure 1: Froh distribution\n")

# --- Supplementary Figure 2: FST distribution ---
supp_fig2 <- ggplot(combined_data, aes(x = FST)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7, color = "black", size = 0.2) +
  geom_vline(xintercept = threshold_995, color = "red", linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = threshold_99, color = "orange", linetype = "dashed", size = 0.6) +
  annotate("text", x = threshold_995, y = Inf, 
           label = sprintf("99.5%% (%.4f)", threshold_995), 
           vjust = 1.5, hjust = -0.1, color = "red", size = 3) +
  annotate("text", x = threshold_99, y = Inf, 
           label = sprintf("99%% (%.4f)", threshold_99), 
           vjust = 3, hjust = -0.1, color = "orange", size = 3) +
  labs(
    title = "Distribution of FST Values",
    x = expression(F[ST]),
    y = "Count"
  ) +
  theme_classic(base_size = base_text_size) +
  theme(
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )

ggsave(paste0(output_prefix, "_supplementary_fst_distribution.png"), 
       supp_fig2, width = 10, height = 6, dpi = figure_dpi)
cat("  ✓ Supplementary Figure 2: FST distribution\n")

# --- Supplementary Figure 3: SNP category breakdown ---
category_counts <- manhattan_data %>%
  count(point_type) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = sprintf("%s\n(n=%d, %.1f%%)", point_type, n, percentage)
  )

supp_fig3 <- ggplot(category_counts, aes(x = "", y = n, fill = point_type)) +
  geom_col(color = "black", size = 0.3) +
  coord_polar("y", start = 0) +
  scale_fill_manual(
    values = c(
      "Neutral" = manhattan_colors$neutral,
      "Private allele only" = manhattan_colors$private_only,
      "FST outlier only" = manhattan_colors$fst_only,
      "FST outlier + Private" = manhattan_colors$both
    ),
    name = "SNP Category"
  ) +
  labs(title = "SNP Classification Summary") +
  theme_void(base_size = base_text_size) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = title_text_size),
    legend.position = "right",
    legend.text = element_text(size = axis_text_size)
  )

ggsave(paste0(output_prefix, "_supplementary_snp_categories.png"), 
       supp_fig3, width = 8, height = 6, dpi = figure_dpi)
cat("  ✓ Supplementary Figure 3: SNP categories pie chart\n")

# --- Supplementary Figure 4: Private alleles per sample ---
supp_fig4 <- ggplot(private_data, 
                    aes(x = reorder(population, -private_alleles_per_sample), 
                        y = private_alleles_per_sample,
                        fill = population)) +
  geom_col(color = "black", size = 0.3) +
  labs(
    title = "Private Alleles per Sample by Population",
    x = "Population",
    y = "Private Alleles per Sample"
  ) +
  scale_fill_viridis_d(option = private_allele_colors) +
  theme_classic(base_size = base_text_size) +
  theme(
    axis.text.x = element_text(angle = axis_text_angle, hjust = 1, size = axis_text_size),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )

ggsave(paste0(output_prefix, "_supplementary_private_per_sample.png"), 
       supp_fig4, width = 10, height = 6, dpi = figure_dpi)
cat("  ✓ Supplementary Figure 4: Private alleles per sample\n\n")

# =============================================================================
# SECTION 9: SUMMARY STATISTICS
# =============================================================================

cat("=============================================================================\n")
cat("FIGURE GENERATION SUMMARY\n")
cat("=============================================================================\n\n")

cat("Main figures created:\n")
cat(sprintf("  ★ %s\n", output_file))
cat(sprintf("  ★ %s (high-res)\n", output_file_hires))
cat("\nIndividual panels:\n")
cat(sprintf("  - %s_panel_A_private_alleles.png\n", output_prefix))
cat(sprintf("  - %s_panel_B_froh.png\n", output_prefix))
cat(sprintf("  - %s_panel_C_manhattan.png\n", output_prefix))
cat("\nSupplementary figures:\n")
cat(sprintf("  - %s_supplementary_froh_distribution.png\n", output_prefix))
cat(sprintf("  - %s_supplementary_fst_distribution.png\n", output_prefix))
cat(sprintf("  - %s_supplementary_snp_categories.png\n", output_prefix))
cat(sprintf("  - %s_supplementary_private_per_sample.png\n", output_prefix))

cat("\n=============================================================================\n")
cat("STATISTICS SUMMARY\n")
cat("=============================================================================\n\n")

# Froh statistics
cat("FROH STATISTICS:\n")
froh_summary <- froh_data %>%
  summarise(
    n_individuals = n(),
    mean_froh = mean(Froh),
    median_froh = median(Froh),
    sd_froh = sd(Froh),
    min_froh = min(Froh),
    max_froh = max(Froh)
  )
print(froh_summary)
cat("\nTop 3 populations by mean Froh:\n")
top_froh <- froh_data %>%
  group_by(population) %>%
  summarise(mean_froh = mean(Froh), .groups = "drop") %>%
  arrange(desc(mean_froh)) %>%
  head(3)
print(top_froh)

# Private allele statistics
cat("\n\nPRIVATE ALLELE STATISTICS:\n")
private_summary <- private_data %>%
  summarise(
    n_populations = n(),
    total_private = sum(private_alleles),
    mean_private = mean(private_alleles),
    median_private = median(private_alleles),
    sd_private = sd(private_alleles)
  )
print(private_summary)
cat("\nTop 3 populations by private alleles:\n")
top_private <- private_data %>%
  arrange(desc(private_alleles)) %>%
  head(3) %>%
  select(population, private_alleles, private_alleles_per_sample)
print(top_private)

# FST statistics
cat("\n\nFST AND OUTLIER STATISTICS:\n")
fst_summary <- combined_data %>%
  summarise(
    n_snps = n(),
    mean_fst = mean(FST, na.rm = TRUE),
    median_fst = median(FST, na.rm = TRUE),
    n_outliers_995 = sum(Outlier_995, na.rm = TRUE),
    n_outliers_99 = sum(Outlier_99, na.rm = TRUE),
    n_private = sum(Is_Private_Any, na.rm = TRUE),
    n_both = sum(Outlier_995 & Is_Private_Any, na.rm = TRUE)
  )
print(fst_summary)

cat("\n\nSNP CATEGORY BREAKDOWN:\n")
print(category_counts %>% select(point_type, n, percentage))

cat("\n=============================================================================\n")
cat("PLOTTING COMPLETE!\n")
cat("=============================================================================\n")