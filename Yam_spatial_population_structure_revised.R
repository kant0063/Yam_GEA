# ===============================================================================
# YAM SPATIAL POPULATION STRUCTURE ANALYSIS WITH TESS3
# Revised version with consistent colors and ordering across all visualizations
# ===============================================================================

# Load required libraries
library(vcfR)
library(tess3r)
library(maps)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(rworldmap)
library(scatterpie)  # For pie chart maps

# ===== LOAD AND PROCESS VCF FILE =====

# Read VCF file
#vcf <- read.vcfR("Calling_ALL_Rotundata_Allc05.vcf.gz")
vcf <- read.vcfR("Yam_no_missing.vcf")

cat("VCF loaded with", nrow(vcf@fix), "SNPs and", ncol(vcf@gt) - 1, "samples\n")

# Extract genotypes as numeric matrix (0, 1, 2 format for diploid)
gt_numeric <- extract.gt(vcf, element = "GT", as.numeric = TRUE)

# Transpose so samples are rows and SNPs are columns (required for TESS3)
genotype <- t(gt_numeric)

# Get sample names
sample_names <- rownames(genotype)

cat("Genotype matrix dimensions:", dim(genotype), "\n")
cat("Samples (rows):", nrow(genotype), "\n")
cat("SNPs (columns):", ncol(genotype), "\n")

# ===== LOAD COORDINATES =====

# Load coordinate data
coordinates <- read.csv("Yam_meta_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Keep necessary columns including population assignment
# Adjust column names based on your actual metadata file
coordinates <- coordinates %>%
  select(SampleID, Longitude, Latitude, Population_assignment) %>%  # Include population column
  rename(Accession = SampleID)

cat("\nCoordinate data loaded with", nrow(coordinates), "accessions\n")
cat("Population assignments found:", length(unique(coordinates$Population_assignment)), "unique populations\n")

# ===== MATCH GENOTYPE AND COORDINATE DATA =====

# Match sample names between VCF and coordinates
common_samples <- intersect(sample_names, coordinates$Accession)
cat("\nCommon samples between VCF and coordinates:", length(common_samples), "\n")

# Subset and order both datasets to match
genotype_matched <- genotype[rownames(genotype) %in% common_samples, ]
coordinates_matched <- coordinates[coordinates$Accession %in% common_samples, ]

# Ensure same order
coordinates_matched <- coordinates_matched[match(rownames(genotype_matched), 
                                                 coordinates_matched$Accession), ]

# Verify alignment
all(rownames(genotype_matched) == coordinates_matched$Accession)

# ===== REMOVE SAMPLES WITH MISSING COORDINATES =====

# Identify samples with missing lat/long
missing_coords <- which(!complete.cases(coordinates_matched[, c("Latitude", "Longitude")]))

if (length(missing_coords) > 0) {
  cat("\nRemoving", length(missing_coords), "samples with missing coordinates\n")
  cat("Sample indices:", missing_coords, "\n")
  
  genotype_final <- genotype_matched[-missing_coords, ]
  coordinates_final <- coordinates_matched[-missing_coords, ]
} else {
  cat("\nNo missing coordinates found\n")
  genotype_final <- genotype_matched
  coordinates_final <- coordinates_matched
}

# ===== HANDLE MISSING GENOTYPE DATA =====

# Impute missing genotypes with column means (per SNP)
genotype_final <- apply(genotype_final, 2, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
})

# Check for any remaining NAs
cat("\nRemaining NAs in genotype matrix:", sum(is.na(genotype_final)), "\n")

# Extract just the coordinate matrix for TESS3
coord_matrix <- as.matrix(coordinates_final[, c("Longitude", "Latitude")])

cat("\nFinal dataset dimensions:\n")
cat("Samples:", nrow(genotype_final), "\n")
cat("SNPs:", ncol(genotype_final), "\n")
cat("Coordinates:", nrow(coord_matrix), "\n")

# ===== PLOT SAMPLING LOCATIONS =====

plot(coord_matrix, pch = 19, cex = 0.5,
     xlab = "Longitude (°E)", ylab = "Latitude (°N)",
     main = "Sampling Locations")
maps::map(add = TRUE, interior = TRUE)

# ===== RUN TESS3 ANALYSIS =====

cat("\nRunning TESS3 analysis for K = 1 to 7...\n")
cat("This may take several minutes...\n")

tess3.obj <- tess3(X = genotype_final, 
                   coord = coord_matrix, 
                   K = 1:7,
                   method = "projected.ls", 
                   ploidy = 2, 
                   openMP.core.num = 4)

# ===== PLOT CROSS-VALIDATION SCORES =====

plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score",
     main = "TESS3 Cross-Validation Scores")

# Extract cross-validation scores
cat("\nExtracting cross-validation scores...\n")

if (!is.null(tess3.obj$crossvalid)) {
  cv_scores <- tess3.obj$crossvalid
} else if (!is.null(tess3.obj$crossentropy)) {
  cv_scores <- tess3.obj$crossentropy
} else {
  cv_scores <- sapply(1:7, function(k) {
    tess3.obj[[k]]$crossentropy
  })
}

cat("\nCross-validation scores:\n")
for (k in 1:length(cv_scores)) {
  cat("K =", k, ": ", cv_scores[k], "\n")
}

# Find optimal K
if (length(cv_scores) > 1) {
  optimal_K <- which.min(cv_scores[2:length(cv_scores)]) + 1
  cat("\nOptimal K (minimum CV score):", optimal_K, "\n")
  cat("CV score at optimal K:", cv_scores[optimal_K], "\n")
} else {
  cat("\nCould not determine optimal K automatically.\n")
  optimal_K <- NA
}

# ===============================================================================
# DEFINE CONSISTENT COLOR PALETTES FOR EACH K VALUE
# ===============================================================================

# Define comprehensive color palettes for each K
# Using distinct, colorblind-friendly colors
color_palettes <- list(
  K2 = c("#E41A1C", "#377EB8"),
  K3 = c("#E41A1C", "#377EB8", "#4DAF4A"),
  K4 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
  K5 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
  K6 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"),
  K7 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628")
)

# Create tess3r-compatible palettes for each K
my.palettes <- list()
for (K in 2:7) {
  my.palettes[[paste0("K", K)]] <- CreatePalette(color_palettes[[paste0("K", K)]], K)
}

# ===============================================================================
# CREATE STRUCTURE-LIKE BARPLOTS WITH CONSISTENT COLORS
# ===============================================================================

# Create ggplot2 barplots for K=2 to K=7
library(gridExtra)

barplot_list <- list()

for (K in 2:7) {
  q.matrix <- qmatrix(tess3.obj, K = K)
  
  # Convert to data frame for ggplot
  q_df_plot <- as.data.frame(q.matrix)
  colnames(q_df_plot) <- paste0("Pop", 1:K)
  q_df_plot$Sample <- 1:nrow(q_df_plot)
  
  # Convert to long format
  q_long_plot <- q_df_plot %>%
    pivot_longer(cols = starts_with("Pop"),
                 names_to = "Population",
                 values_to = "Proportion")
  
  # Create barplot
  p <- ggplot(q_long_plot, aes(x = Sample, y = Proportion, fill = Population)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    scale_fill_manual(values = color_palettes[[paste0("K", K)]]) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = paste("Ancestry matrix (K =", K, ")"),
         x = "Individuals",
         y = "Ancestry proportions") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  barplot_list[[paste0("K", K)]] <- p
}

# Display in grid
grid.arrange(grobs = barplot_list, ncol = 2, nrow = 3)

# ===== SAVE COMBINED BARPLOT =====

# Recreate barplot list for saving
barplot_list <- list()

for (K in 2:7) {
  q.matrix <- qmatrix(tess3.obj, K = K)
  
  # Convert to data frame for ggplot
  q_df_plot <- as.data.frame(q.matrix)
  colnames(q_df_plot) <- paste0("Pop", 1:K)
  q_df_plot$Sample <- 1:nrow(q_df_plot)
  
  # Convert to long format
  q_long_plot <- q_df_plot %>%
    pivot_longer(cols = starts_with("Pop"),
                 names_to = "Population",
                 values_to = "Proportion")
  
  # Create barplot
  p <- ggplot(q_long_plot, aes(x = Sample, y = Proportion, fill = Population)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    scale_fill_manual(values = color_palettes[[paste0("K", K)]]) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = paste("Ancestry matrix (K =", K, ")"),
         x = "Individuals",
         y = "Ancestry proportions") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  barplot_list[[paste0("K", K)]] <- p
}

# Combine and save
combined_barplots <- wrap_plots(barplot_list, ncol = 2, nrow = 3)

ggsave("tess3_barplots_K2_to_K7_consistent.png", 
       combined_barplots,
       width = 16, height = 20, dpi = 300)

cat("Combined barplot saved to: tess3_barplots_K2_to_K7_consistent.png\n")

# ===== SAVE INDIVIDUAL FILES WITH CONSISTENT COLORS =====

for (K in 2:7) {
  q.matrix <- qmatrix(tess3.obj, K = K)
  
  # Convert to data frame for ggplot
  q_df_plot <- as.data.frame(q.matrix)
  colnames(q_df_plot) <- paste0("Pop", 1:K)
  q_df_plot$Sample <- 1:nrow(q_df_plot)
  
  # Convert to long format
  q_long_plot <- q_df_plot %>%
    pivot_longer(cols = starts_with("Pop"),
                 names_to = "Population",
                 values_to = "Proportion")
  
  # Create barplot
  p <- ggplot(q_long_plot, aes(x = Sample, y = Proportion, fill = Population)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    scale_fill_manual(values = color_palettes[[paste0("K", K)]]) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = paste("Ancestry matrix (K =", K, ")"),
         x = "Individuals",
         y = "Ancestry proportions") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  # Save
  ggsave(paste0("tess3_barplot_K", K, "_consistent.png"),
         p,
         width = 14, height = 7, dpi = 300)
  
  cat("Saved: tess3_barplot_K", K, "_consistent.png\n")
}

# ===============================================================================
# CREATE SPATIAL MAPS WITH CONSISTENT COLORS
# ===============================================================================

# Base R maps with consistent colors
png("tess3_spatial_maps_baseR_consistent.png", width = 16, height = 20, units = "in", res = 300)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))

for (K in 2:7) {
  q.matrix <- qmatrix(tess3.obj, K = K)
  
  plot(q.matrix, coord_matrix, method = "map.max", 
       interpol = FieldsKrigModel(10),
       main = paste("Ancestry coefficients (K =", K, ")"),
       xlab = "Longitude", ylab = "Latitude",
       resolution = c(300, 300), cex = 0.4,
       col.palette = my.palettes[[paste0("K", K)]],
       cex.main = 1.3,
       cex.lab = 1.2,
       cex.axis = 1.1)
}

par(mfrow = c(1, 1))
dev.off()

cat("Base R spatial maps saved to 'tess3_spatial_maps_baseR_consistent.png'\n")

# ===============================================================================
# CREATE GGPLOT MAPS WITH CONSISTENT COLORS
# ===============================================================================

map.polygon <- getMap(resolution = "low")

plot_list <- list()

for (K in 2:7) {
  q.matrix <- qmatrix(tess3.obj, K = K)
  
  pl <- ggtess3Q(q.matrix, coord_matrix, map.polygon = map.polygon,
                 col.palette = my.palettes[[paste0("K", K)]])
  
  p <- pl +
    geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
    xlim(-20, 15) +
    ylim(-10, 15) +
    coord_equal() +
    geom_point(data = data.frame(Longitude = coord_matrix[, 1], 
                                 Latitude = coord_matrix[, 2]),
               aes(x = Longitude, y = Latitude),
               size = 0.5, color = "black") +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_bw() +
    ggtitle(paste("Ancestry coefficients (K =", K, ")")) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11),
      legend.key.size = unit(0.6, "cm")
    )
  
  plot_list[[paste0("K", K)]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 2, nrow = 3)

ggsave("Yam_tess3_spatial_maps_ggplot_consistent.png", combined_plot, 
       width = 18, height = 22, dpi = 300)

cat("ggplot spatial maps saved to 'Yam_tess3_spatial_maps_ggplot_consistent.png'\n")

print(combined_plot)

# ===============================================================================
# CREATE PIE CHART MAP FOR OPTIMAL K
# ===============================================================================

# Use optimal K or default to K=5 if not determined
K_for_pie <- ifelse(is.na(optimal_K), 5, optimal_K)

cat("\nCreating pie chart map for K =", K_for_pie, "\n")

# Extract Q-matrix for chosen K
q.matrix <- qmatrix(tess3.obj, K = K_for_pie)

# Prepare data for scatterpie
q_df <- as.data.frame(q.matrix, stringsAsFactors = FALSE)

# Get actual population names from metadata
# Match samples to get their population assignments
sample_pop_map <- data.frame(
  Sample = rownames(genotype_final),
  stringsAsFactors = FALSE
) %>%
  left_join(coordinates_final %>% 
              mutate(Sample = Accession) %>%
              select(Sample, Population_assignment),
            by = "Sample")

# Get unique populations and create a mapping to K clusters
# Assign cluster names based on which population is most common in each cluster
pop_assignments <- character(K_for_pie)
for (k in 1:K_for_pie) {
  # Find samples with highest membership in this cluster
  top_samples <- order(q.matrix[, k], decreasing = TRUE)[1:min(10, nrow(q.matrix))]
  # Get most common population among these samples
  common_pop <- names(sort(table(sample_pop_map$Population_assignment[top_samples]), 
                           decreasing = TRUE))[1]
  pop_assignments[k] <- common_pop
}

# Handle duplicates by adding numbers
pop_assignments_unique <- make.unique(pop_assignments, sep = "_")

colnames(q_df) <- pop_assignments_unique

# Add Sample as a proper column using rownames from genotype_final
if (is.null(rownames(q.matrix)) || length(rownames(q.matrix)) == 0) {
  sample_names_for_df <- rownames(genotype_final)
} else {
  sample_names_for_df <- rownames(q.matrix)
}

q_df$Sample <- sample_names_for_df
q_df$Longitude <- coord_matrix[, 1]
q_df$Latitude <- coord_matrix[, 2]

# Add actual population assignment from metadata
q_df <- q_df %>%
  left_join(sample_pop_map, by = "Sample")

# Verify Sample column exists
cat("Columns in q_df:", paste(colnames(q_df), collapse = ", "), "\n")
cat("Number of rows in q_df:", nrow(q_df), "\n")
if (nrow(q_df) > 0) {
  cat("First few Sample values:", paste(head(q_df$Sample), collapse = ", "), "\n")
}

# Get world map data
world_map <- map_data("world")

# Define population order and colors for this K
pop_names <- pop_assignments_unique
pop_colors <- setNames(color_palettes[[paste0("K", K_for_pie)]], pop_names)

# Create pie chart map
pie_map <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray60", linewidth = 0.3) +
  geom_scatterpie(aes(x = Longitude, y = Latitude),
                  data = q_df,
                  cols = pop_names,
                  color = NA, 
                  alpha = 0.8,
                  pie_scale = 1.5) +
  scale_fill_manual(values = pop_colors, name = "Population") +
  coord_fixed(xlim = c(min(coord_matrix[, 1]) - 5, max(coord_matrix[, 1]) + 5),
              ylim = c(min(coord_matrix[, 2]) - 5, max(coord_matrix[, 2]) + 5)) +
  labs(title = paste0("Population Structure (K = ", K_for_pie, ")"),
       x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# Display the map
print(pie_map)

# Save the map with larger dimensions
ggsave(paste0("Yam_pie_chart_map_K", K_for_pie, ".png"), 
       pie_map, 
       width = 16, height = 10, dpi = 300)

cat("Pie chart map saved to: Yam_pie_chart_map_K", K_for_pie, ".png\n")

# ===============================================================================
# CREATE BARPLOTS ORGANIZED BY GEOGRAPHIC REGION
# ===============================================================================

# Use actual region/population data from metadata
# Replace the longitude-based regions with actual Population_assignment from metadata
q_df$Region <- q_df$Population_assignment

# Determine dominant population for each sample (for sorting)
pop_cols <- pop_assignments_unique
q_df$DominantPop <- pop_cols[apply(q_df[, pop_cols], 1, which.max)]

# Convert to long format for ggplot
q_long <- q_df %>%
  select(Sample, Longitude, Latitude, Region, DominantPop, all_of(pop_cols)) %>%
  pivot_longer(cols = all_of(pop_cols),
               names_to = "Population",
               values_to = "Proportion")

# Convert Population to factor with specific order (for consistent colors)
q_long$Population <- factor(q_long$Population, levels = pop_names)

# Create sample order: first by region, then by dominant population
sample_order <- q_df %>%
  arrange(Region, DominantPop) %>%
  pull(Sample)

# Apply ordering
q_long$Sample <- factor(q_long$Sample, levels = sample_order)

# Create stacked barplot organized by region (actual populations)
# Use facet_wrap with scales and space to handle long labels better
region_barplot <- ggplot(q_long, aes(x = Sample, y = Proportion, fill = Population)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  scale_fill_manual(values = pop_colors, 
                    name = "Ancestry Cluster",
                    labels = paste0("Cluster ", 1:K_for_pie, "\n(", pop_names, ")")) +
  facet_grid(~ Region, scales = "free_x", space = "free_x") +
  labs(title = paste0("Population Structure by Population Assignment (K = ", K_for_pie, ")"),
       x = "Individuals",
       y = "Ancestry Proportion") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm"),
    legend.position = "right",
    strip.text = element_text(size = 11, face = "bold", color = "black", 
                              margin = margin(t = 5, b = 5)),
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.5),
    panel.spacing = unit(0.3, "lines")
  )

# Display the plot
print(region_barplot)

# Save the plot with larger dimensions to ensure facet labels are readable
ggsave(paste0("Yam_barplot_by_population_K", K_for_pie, ".png"),
       region_barplot,
       width = 24, height = 8, dpi = 300)

cat("Population barplot saved to: Yam_barplot_by_population_K", K_for_pie, ".png\n")

# ===== CREATE ALTERNATIVE VERSION WITH ABBREVIATED LABELS =====

# Create abbreviated population names for facets (for readability)
pop_abbrev <- data.frame(
  Region = unique(q_df$Region),
  RegionAbbrev = make.names(abbreviate(unique(q_df$Region), minlength = 8))
)

q_long_abbrev <- q_long %>%
  left_join(pop_abbrev, by = "Region")

# Reorder by abbreviated names
sample_order_abbrev <- q_df %>%
  left_join(pop_abbrev, by = "Region") %>%
  arrange(RegionAbbrev, DominantPop) %>%
  pull(Sample)

q_long_abbrev$Sample <- factor(q_long_abbrev$Sample, levels = sample_order_abbrev)

# Create barplot with abbreviated facet labels
region_barplot_abbrev <- ggplot(q_long_abbrev, aes(x = Sample, y = Proportion, fill = Population)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  scale_fill_manual(values = pop_colors, 
                    name = "Ancestry Cluster",
                    labels = paste0("K", 1:K_for_pie, ": ", pop_names)) +
  facet_grid(~ RegionAbbrev, scales = "free_x", space = "free_x") +
  labs(title = paste0("Population Structure by Population Assignment (K = ", K_for_pie, ")"),
       subtitle = "Facet labels show abbreviated population names",
       x = "Individuals",
       y = "Ancestry Proportion") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "right",
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.5),
    panel.spacing = unit(0.3, "lines")
  )

# Save abbreviated version
ggsave(paste0("Yam_barplot_by_population_K", K_for_pie, "_abbreviated.png"),
       region_barplot_abbrev,
       width = 24, height = 8, dpi = 300)

cat("Abbreviated barplot saved to: Yam_barplot_by_population_K", K_for_pie, "_abbreviated.png\n")

# ===============================================================================
# CREATE COMBINED PIE CHART AND INTERPOLATION MAP FOR K=7
# ===============================================================================

cat("\nCreating combined pie chart and interpolation map for K =", K_for_pie, "\n")

# Create interpolation map for the optimal K using ggplot
q.matrix_for_interp <- qmatrix(tess3.obj, K = K_for_pie)

# Create the interpolation map using ggtess3Q
interp_map <- ggtess3Q(q.matrix_for_interp, coord_matrix, map.polygon = map.polygon,
                       col.palette = my.palettes[[paste0("K", K_for_pie)]])

interp_map <- interp_map +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  xlim(min(coord_matrix[, 1]) - 5, max(coord_matrix[, 1]) + 5) +
  ylim(min(coord_matrix[, 2]) - 5, max(coord_matrix[, 2]) + 5) +
  coord_equal() +
  geom_point(data = data.frame(Longitude = coord_matrix[, 1], 
                               Latitude = coord_matrix[, 2]),
             aes(x = Longitude, y = Latitude),
             size = 1, color = "black") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  ggtitle(paste0("Interpolated Ancestry Coefficients (K = ", K_for_pie, ")")) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.6, "cm"),
    legend.position = "right"
  )

# Recreate pie chart map with matching theme
pie_map_combined <- ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray60", linewidth = 0.3) +
  geom_scatterpie(aes(x = Longitude, y = Latitude),
                  data = q_df,
                  cols = pop_names,
                  color = NA, 
                  alpha = 0.8,
                  pie_scale = 1.5) +
  scale_fill_manual(values = pop_colors, 
                    name = "Ancestry Cluster",
                    labels = paste0("K", 1:K_for_pie, ": ", pop_names)) +
  coord_fixed(xlim = c(min(coord_matrix[, 1]) - 5, max(coord_matrix[, 1]) + 5),
              ylim = c(min(coord_matrix[, 2]) - 5, max(coord_matrix[, 2]) + 5)) +
  labs(title = paste0("Individual Ancestry Proportions (K = ", K_for_pie, ")"),
       x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.6, "cm"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# Combine the two maps side by side using patchwork
combined_maps <- pie_map_combined + interp_map +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = paste0("Population Structure Visualization (K = ", K_for_pie, ")"),
    subtitle = "Left: Individual pie charts | Right: Interpolated ancestry surface",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30")
    )
  )

# Display the combined plot
print(combined_maps)

# Save the combined map
ggsave(paste0("Yam_combined_pie_and_interpolation_K", K_for_pie, ".png"),
       combined_maps,
       width = 20, height = 10, dpi = 300)

cat("Combined pie chart and interpolation map saved to: Yam_combined_pie_and_interpolation_K", K_for_pie, ".png\n")

# ===============================================================================
# CREATE COMBINED FIGURE: MAP + BARPLOT
# ===============================================================================

combined_map_bar <- pie_map / region_barplot + 
  plot_layout(heights = c(2, 1))

# Save combined figure with larger size
ggsave(paste0("Yam_map_and_barplot_K", K_for_pie, ".png"),
       combined_map_bar,
       width = 20, height = 16, dpi = 300)

cat("Combined map and barplot saved to: Yam_map_and_barplot_K", K_for_pie, ".png\n")

# ===============================================================================
# SAVE Q-MATRICES FOR FURTHER ANALYSIS
# ===============================================================================

for (K in 2:7) {
  q.matrix <- qmatrix(tess3.obj, K = K)
  q_df_export <- data.frame(
    Sample = rownames(genotype_final),
    q.matrix
  )
  colnames(q_df_export) <- c("Sample", paste0("Pop", 1:K))
  
  # Add coordinates
  q_df_export <- q_df_export %>%
    left_join(coordinates_final %>% 
                mutate(Sample = Accession) %>%
                select(Sample, Longitude, Latitude),
              by = "Sample")
  
  write.csv(q_df_export, paste0("tess3_qmatrix_K", K, "_with_coords.csv"), row.names = FALSE)
}

cat("\nQ-matrices with coordinates saved for K = 2 to 7\n")

# ===============================================================================
# FIGURE 1 PANELS B AND C (K = 7) FOR THE COMPOSITE BUILT IN SCRIPT 17
# ===============================================================================
# Panels are drawn without titles (panel letters are added in the composite) at K = 7,
# the K reported in the manuscript, whatever the cross-validation optimum. They are saved
# as ggplot objects (results/strict/Fig1BC_panels.rds), which script 17 combines with the
# PCA (panel A), and as a stand-alone PNG (figures/strict/Fig1BC_tess3.png).
K_fig <- 7
q_fig <- qmatrix(tess3.obj, K = K_fig)
fig_names <- character(K_fig)
for (k in 1:K_fig) {                      # name each cluster after the group most common among its top-10 members
  top <- order(q_fig[, k], decreasing = TRUE)[1:min(10, nrow(q_fig))]
  fig_names[k] <- names(sort(table(sample_pop_map$Population_assignment[top]), decreasing = TRUE))[1]
}
fig_names <- make.unique(fig_names, sep = " ")
df_fig <- as.data.frame(q_fig); colnames(df_fig) <- fig_names
df_fig$Longitude <- coord_matrix[, 1]; df_fig$Latitude <- coord_matrix[, 2]
fig_cols <- setNames(color_palettes[[paste0("K", K_fig)]], fig_names)
xl_fig <- range(coord_matrix[, 1]) + c(-1.5, 1.5); yl_fig <- range(coord_matrix[, 2]) + c(-1.5, 1.5)

fig1B <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray92", color = "gray55", linewidth = 0.25) +
  geom_scatterpie(aes(x = Longitude, y = Latitude), data = df_fig, cols = fig_names, color = NA, alpha = 0.9, pie_scale = 0.55) +
  scale_fill_manual(values = fig_cols, name = paste0("Ancestry cluster (K = ", K_fig, ")")) +
  coord_fixed(xlim = xl_fig, ylim = yl_fig, expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") + theme_bw() +
  theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8), panel.grid.minor = element_blank())

fig1C <- ggtess3Q(q_fig, coord_matrix, map.polygon = map.polygon, col.palette = my.palettes[[paste0("K", K_fig)]]) +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group), linewidth = 0.25, colour = "gray30") +
  geom_point(data = data.frame(Longitude = coord_matrix[, 1], Latitude = coord_matrix[, 2]),
             aes(x = Longitude, y = Latitude), size = 0.6, colour = "black") +
  coord_equal(xlim = xl_fig, ylim = yl_fig, expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") + theme_bw() + theme(panel.grid.minor = element_blank())

dir.create(file.path("results", "strict"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("figures", "strict"), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(B = fig1B, C = fig1C, K = K_fig, cluster_names = fig_names, colours = fig_cols),
        file.path("results", "strict", "Fig1BC_panels.rds"))
ggsave(file.path("figures", "strict", "Fig1BC_tess3.png"),
       fig1B + fig1C + plot_layout(ncol = 2, guides = "collect"), width = 12, height = 4.5, dpi = 300, bg = "white")
cat("Figure 1 panels B-C saved to results/strict/Fig1BC_panels.rds and figures/strict/Fig1BC_tess3.png\n")
cat("Cluster names (K = 7):", paste(fig_names, collapse = "; "), "\n")

# ===============================================================================
# DATA FOR FIGURES S2 (ancestry bar plots, K = 2-7) AND S3 (interpolated maps, K = 2-7)
# ===============================================================================
# Script 17 draws the composite Figures S2 and S3 from this file (panels A-F = K 2-7),
# with the same colours as above and accessions grouped by genetic group.
tess3_store <- list(Q = setNames(lapply(2:7, function(K) qmatrix(tess3.obj, K = K)), paste0("K", 2:7)),
                    sample = sample_pop_map$Sample, group = sample_pop_map$Population_assignment,
                    coords = coord_matrix, bar_palettes = color_palettes[paste0("K", 2:7)],
                    map_palettes = my.palettes[paste0("K", 2:7)], map.polygon = map.polygon)
saveRDS(tess3_store, file.path("results", "strict", "tess3_K2_7.rds"))
cat("Figure S2/S3 data saved to results/strict/tess3_K2_7.rds\n")

# ===============================================================================
# SUMMARY OUTPUT
# ===============================================================================

cat("\n=== TESS3 ANALYSIS COMPLETE ===\n")
cat("Final sample size:", nrow(genotype_final), "\n")
cat("Number of SNPs:", ncol(genotype_final), "\n")
cat("Optimal K:", optimal_K, "\n")
cat("\nFiles created with CONSISTENT COLORS:\n")
cat("- tess3_barplots_K2_to_K7_consistent.png\n")
cat("- tess3_barplot_K*_consistent.png (individual plots)\n")
cat("- tess3_spatial_maps_baseR_consistent.png\n")
cat("- Yam_tess3_spatial_maps_ggplot_consistent.png\n")
cat("- Yam_pie_chart_map_K*.png\n")
cat("- Yam_barplot_by_region_K*.png\n")
cat("- Yam_map_and_barplot_K*.png (combined figure)\n")
cat("- tess3_qmatrix_K*_with_coords.csv\n")
cat("- results/strict/Fig1BC_panels.rds + figures/strict/Fig1BC_tess3.png (Figure 1 B-C; composite made by script 17)\n")

cat("\n=== COLOR CONSISTENCY FEATURES ===\n")
cat("✓ Same colors used across all K values\n")
cat("✓ Consistent population order in all plots\n")
cat("✓ Color palettes defined at beginning of script\n")
cat("✓ Samples sorted by dominant population in barplots\n")
cat("✓ Regional organization available in barplots\n")