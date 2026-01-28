# ==========================================================
# Combined RDA and Genomic Offset Two-Panel Figure
# Panel A: RDA biplot
# Panel B: Regional genomic offset map
# FIXED VERSION: Improved sample matching and debugging
# ==========================================================

# ---------------------- LIBRARIES ----------------------
required_pkgs <- c("vcfR", "vegan", "ggplot2", "dplyr", "LEA", "lfmm", 
                   "terra", "geodata", "fields", "maps", "cowplot", "gridExtra")
missing <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(missing)) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing)
}

library(vcfR)
library(vegan)
library(ggplot2)
library(dplyr)
library(LEA)
library(lfmm)
library(terra)
library(geodata)
library(fields)
library(maps)
library(cowplot)
library(gridExtra)

# ---------------------- USER PARAMETERS ----------------------
# VCF and environmental data
vcf_file <- "Yam_no_missing.vcf"
env_csv <- "Yam_bioclim_data_extracted_V3.csv"
meta_file <- "Yam_meta_data.csv"

# Genomic offset parameters
res_minutes <- 10
cmip6_models <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5",
                  "CNRM-CM6-1", "IPSL-CM6A-LR", "MIROC6")
ssp_scenarios <- c("245", "585")
future_periods <- c("2041-2060", "2081-2100")
nc <- 200
K_lfmm <- 7
pvalue_cutoff <- 1e-5
work_dir <- file.path(getwd(), "cmip6_cache")
out_dir <- file.path(getwd(), "combined_figure_output")

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------- PART 1: RDA ANALYSIS ----------------------
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART 1: RDA ANALYSIS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Read VCF and extract genotypes
cat("1. Reading VCF file...\n")
vcf <- read.vcfR(vcf_file)
gt_numeric <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
cat("   VCF samples:", ncol(gt_numeric), "\n")
cat("   VCF SNPs:", nrow(gt_numeric), "\n")

# Load environmental data
cat("\n2. Loading environmental data...\n")
env_data <- read.csv(env_csv, header = TRUE, stringsAsFactors = FALSE)
cat("   Environmental data rows:", nrow(env_data), "\n")
cat("   Environmental data columns:", ncol(env_data), "\n")

# ---------------------- CRITICAL: SAMPLE MATCHING ----------------------
cat("\n3. Matching samples between VCF and environmental data...\n")

# Get sample names from both sources
vcf_samples <- colnames(gt_numeric)
env_samples <- env_data$SampleID

cat("   VCF samples:", length(vcf_samples), "\n")
cat("   ENV samples:", length(env_samples), "\n")

# Find common samples
common_samples <- intersect(vcf_samples, env_samples)
cat("   Common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  cat("\n ERROR: No common samples found!\n")
  cat("\n First 5 VCF samples:\n")
  print(head(vcf_samples, 5))
  cat("\n First 5 ENV samples:\n")
  print(head(env_samples, 5))
  stop("Sample names do not match between VCF and environmental data")
}

if (length(common_samples) < 10) {
  warning("Very few samples matched (", length(common_samples), "). Check sample names!")
}

# Subset and order both datasets to match
gt_matched <- gt_numeric[, common_samples]
env_matched <- env_data[match(common_samples, env_data$SampleID), ]

cat("   Matched genotype dims:", dim(gt_matched), "\n")
cat("   Matched env dims:", dim(env_matched), "\n")

# Verify order matches
if (!all(colnames(gt_matched) == env_matched$SampleID)) {
  stop("Sample order mismatch after matching!")
}

# Prepare genotype matrix
cat("\n4. Preparing genotype matrix...\n")
cat("   Input gt_matched dims (SNPs × samples):", dim(gt_matched), "\n")

# Impute missing values - apply works across rows (SNPs)
cat("   Imputing missing genotypes...\n")
missing_before <- sum(is.na(gt_matched))

# For each SNP (row), replace NA with mean across samples
gt_imputed <- apply(gt_matched, 1, function(snp) {
  ifelse(is.na(snp), mean(snp, na.rm = TRUE), snp)
})

cat("   After apply dims (samples × SNPs):", dim(gt_imputed), "\n")

# apply() with MARGIN=1 transposes: input is SNPs×samples, output is samples×SNPs
# This is exactly what we need for RDA!
genotype_data <- gt_imputed

missing_after <- sum(is.na(genotype_data))

cat("   Missing before imputation:", missing_before, "\n")
cat("   Missing after imputation:", missing_after, "\n")
cat("   Final genotype matrix dims (samples × SNPs):", dim(genotype_data), "\n")

# Add rownames for samples
rownames(genotype_data) <- common_samples
cat("   Rownames added for", nrow(genotype_data), "samples\n")

# CRITICAL CHECK: Verify dimensions are correct
if (nrow(genotype_data) != length(common_samples)) {
  cat("\n   ERROR: Genotype matrix has wrong number of samples!\n")
  cat("   Expected samples:", length(common_samples), "\n")
  cat("   Actual rows:", nrow(genotype_data), "\n")
  cat("   This suggests a transpose error.\n")
  stop("Genotype matrix dimension error - check transpose operations")
}
if (ncol(genotype_data) < 1000) {
  warning("Very few SNPs (", ncol(genotype_data), "). Expected more for typical genomic analysis.")
}
cat("   ✓ Genotype matrix dimensions verified\n")

# ---------------------- PREPARE ENVIRONMENTAL VARIABLES ----------------------
cat("\n5. Preparing environmental variables...\n")

# Extract only bioclim variables (exclude SampleID, Lat, Lon)
exclude_cols <- c("SampleID", "Longitude", "Latitude")
env_vars <- env_matched[, !colnames(env_matched) %in% exclude_cols, drop = FALSE]

cat("   Environmental variables available:", ncol(env_vars), "\n")
cat("   Variables:", paste(colnames(env_vars), collapse = ", "), "\n")

# Remove variables with no variation or too many missing values
cat("\n6. Cleaning environmental variables...\n")
env_vars_clean <- env_vars[, sapply(env_vars, function(x) {
  # Keep if: has variance AND has at least 50% non-missing
  has_variance <- var(x, na.rm = TRUE) > 0
  enough_data <- sum(!is.na(x)) > (nrow(env_vars) * 0.5)
  has_variance && enough_data
})]

cat("   Variables after cleaning:", ncol(env_vars_clean), "\n")

if (ncol(env_vars_clean) == 0) {
  stop("No valid environmental variables remaining after cleaning!")
}

# Impute missing values with column means
cat("   Imputing missing environmental values...\n")
env_vars_clean <- as.data.frame(lapply(env_vars_clean, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}))

# Remove highly correlated variables
cat("\n7. Checking for highly correlated variables...\n")
env_cor <- cor(env_vars_clean, use = "complete.obs")

identify_highly_correlated <- function(cor_matrix, cutoff = 0.9) {
  upper_tri <- upper.tri(cor_matrix)
  high_cor_pairs <- which(abs(cor_matrix) > cutoff & upper_tri, arr.ind = TRUE)
  
  if (nrow(high_cor_pairs) == 0) return(integer(0))
  
  vars_to_remove <- c()
  for (i in 1:nrow(high_cor_pairs)) {
    var1 <- high_cor_pairs[i, 1]
    var2 <- high_cor_pairs[i, 2]
    mean_cor_var1 <- mean(abs(cor_matrix[var1, -var1]), na.rm = TRUE)
    mean_cor_var2 <- mean(abs(cor_matrix[var2, -var2]), na.rm = TRUE)
    
    if (mean_cor_var1 > mean_cor_var2) {
      vars_to_remove <- c(vars_to_remove, var1)
    } else {
      vars_to_remove <- c(vars_to_remove, var2)
    }
  }
  return(unique(vars_to_remove))
}

highly_corr <- identify_highly_correlated(env_cor, cutoff = 0.9)

if (length(highly_corr) > 0) {
  cat("   Removing", length(highly_corr), "highly correlated variables:\n")
  cat("   ", paste(colnames(env_vars_clean)[highly_corr], collapse = ", "), "\n")
  env_vars_final <- env_vars_clean[, -highly_corr, drop = FALSE]
} else {
  cat("   No highly correlated variables found (|r| > 0.9)\n")
  env_vars_final <- env_vars_clean
}

cat("   Final environmental variables:", ncol(env_vars_final), "\n")
cat("   ", paste(colnames(env_vars_final), collapse = ", "), "\n")

# Scale environmental variables
cat("\n8. Scaling environmental variables...\n")
env_vars_scaled <- scale(env_vars_final)

# Ensure rownames are preserved
rownames(env_vars_scaled) <- common_samples

# ---------------------- FINAL VERIFICATION ----------------------
cat("\n9. Final verification before RDA...\n")
cat("   Genotype data dimensions (samples × SNPs):", dim(genotype_data), "\n")
cat("   Environmental data dimensions (samples × vars):", dim(env_vars_scaled), "\n")
cat("   Genotype sample names (first 3):", head(rownames(genotype_data), 3), "\n")
cat("   Env sample names (first 3):", head(rownames(env_vars_scaled), 3), "\n")

if (nrow(genotype_data) != nrow(env_vars_scaled)) {
  cat("\n ERROR: Dimension mismatch!\n")
  cat("   Genotype samples:", nrow(genotype_data), "\n")
  cat("   Environment samples:", nrow(env_vars_scaled), "\n")
  stop("Cannot proceed with RDA - sample counts don't match")
}

# ---------------------- RUN RDA ----------------------
cat("\n10. Running RDA analysis...\n")
rda_result <- rda(genotype_data ~ ., data = as.data.frame(env_vars_scaled))

# Extract variance explained
total_var <- rda_result$tot.chi
constrained_var <- rda_result$CCA$tot.chi
var_explained <- constrained_var / total_var

cat("   Total variance:", round(total_var, 2), "\n")
cat("   Constrained variance:", round(constrained_var, 2), "\n")
cat("   Variance explained:", round(var_explained * 100, 2), "%\n")

# Extract RDA scores
site_scores <- scores(rda_result, choices = 1:2, display = "sites")
env_scores <- scores(rda_result, choices = 1:2, display = "bp")

# Create data frames for plotting
rda_df <- data.frame(
  Sample = rownames(site_scores),
  RDA1 = site_scores[, 1],
  RDA2 = site_scores[, 2]
)

env_arrows <- data.frame(
  Variable = rownames(env_scores),
  RDA1 = env_scores[, 1],
  RDA2 = env_scores[, 2]
)

# Calculate variance explained by each axis
axis_var <- eigenvals(rda_result) / total_var * 100

# Load metadata for coloring
cat("\n11. Loading metadata for plotting...\n")
if (file.exists(meta_file)) {
  meta <- read.csv(meta_file)
  # Match metadata to RDA samples
  rda_df$Population <- meta$Population.assignment[match(rda_df$Sample, meta$SampleID)]
  cat("   Loaded population data for", sum(!is.na(rda_df$Population)), "samples\n")
} else {
  cat("   Metadata file not found, using default labels\n")
  rda_df$Population <- "All Samples"
}

# Replace NA populations with "Unknown"
rda_df$Population[is.na(rda_df$Population)] <- "Unknown"

# Create RDA plot (Panel A)
cat("\n12. Creating RDA plot (Panel A)...\n")
panel_A <- ggplot(rda_df, aes(x = RDA1, y = RDA2, color = Population)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_segment(data = env_arrows, 
               aes(x = 0, y = 0, xend = RDA1 * 3, yend = RDA2 * 3),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "red", size = 1, alpha = 0.8,
               inherit.aes = FALSE) +
  geom_text(data = env_arrows, 
            aes(x = RDA1 * 3.2, y = RDA2 * 3.2, label = Variable),
            color = "black", size = 3, fontface = "bold",
            inherit.aes = FALSE) +
  labs(
    title = "A) RDA: Genomic Variation ~ Environment",
    x = paste0("RDA1 (", round(axis_var[1], 2), "%)"),
    y = paste0("RDA2 (", round(axis_var[2], 2), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  coord_fixed()

cat("   Panel A created successfully\n")

# ---------------------- PART 2: GENOMIC OFFSET ANALYSIS ----------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 2: GENOMIC OFFSET ANALYSIS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Utility functions
normalize_bio_names <- function(n) {
  n2 <- tolower(n)
  digits <- gsub(".*?([0-9]+).*", "\\1", n2)
  digits[grepl("^[^0-9]+$", digits)] <- NA
  idx <- suppressWarnings(as.integer(digits))
  idx[is.na(idx)] <- NA
  paste0("bio", idx)
}

to_spatraster <- function(x) {
  if (inherits(x, "SpatRaster")) return(x)
  rast(x)
}

is_valid_map <- function(m) {
  is.matrix(m) && any(is.finite(m))
}

# Read VCF for LFMM (using LEA format)
cat("1. Preparing genotype data for LFMM...\n")
geno_file <- vcf2geno(vcf_file)
genotype_lea <- LEA::read.geno(geno_file)

cat("   Original genotype dims:", dim(genotype_lea), "\n")

# Remove monomorphic sites
is_polymorphic <- function(g) {
  g2 <- g[g != 9]
  if (length(g2) == 0) return(FALSE)
  u <- unique(g2)
  length(u) > 1
}

keep <- apply(genotype_lea, 2, is_polymorphic)
genotype_lea <- genotype_lea[, keep, drop = FALSE]

cat("   Filtered genotype dims:", dim(genotype_lea), "\n")
cat("   Removed", sum(!keep), "monomorphic SNPs\n")

# Read sample coordinates - USE THE MATCHED ENV DATA
cat("\n2. Extracting sample coordinates...\n")
coords <- data.frame(
  Longitude = env_matched$Longitude, 
  Latitude = env_matched$Latitude
)
rownames(coords) <- env_matched$SampleID

cat("   Coordinate range:\n")
cat("   Longitude: [", min(coords$Longitude), ",", max(coords$Longitude), "]\n")
cat("   Latitude: [", min(coords$Latitude), ",", max(coords$Latitude), "]\n")

# Compute data extent
min_lon <- min(coords$Longitude, na.rm = TRUE)
max_lon <- max(coords$Longitude, na.rm = TRUE)
min_lat <- min(coords$Latitude, na.rm = TRUE)
max_lat <- max(coords$Latitude, na.rm = TRUE)

# Download historical WorldClim
cat("\n3. Downloading/reading historical WorldClim...\n")
clim_hist_raw <- worldclim_global(var = "bio", res = res_minutes, path = work_dir, download = TRUE)
clim_hist <- clim_hist_raw
names(clim_hist) <- paste0("bio", seq_len(nlyr(clim_hist)))
cat("   Historical layers:", paste(names(clim_hist), collapse = ", "), "\n")

# Extract training environment
cat("\n4. Extracting environmental data at sample points...\n")
pts <- data.frame(x = coords$Longitude, y = coords$Latitude)
X.env_df <- terra::extract(clim_hist, pts)
X.env <- X.env_df[, -1, drop = FALSE]

# Impute missing values
X.env <- as.data.frame(lapply(X.env, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x 
}))

# Keep only non-zero variance columns
X.env <- X.env[, sapply(X.env, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

if (ncol(X.env) == 0) {
  stop("No non-zero-variance historic env cols.")
}

X.env_mat <- as.matrix(X.env)
X.env_scaled <- scale(X.env_mat)

cat("   Environmental dims (samples × vars):", dim(X.env_mat), "\n")
cat("   Variables:", paste(colnames(X.env), collapse = ", "), "\n")

# Run LFMM
cat("\n5. Running LFMM (this may take a while)...\n")
geno_filtered <- genotype_lea[, colSums(genotype_lea == 9) == 0]
cat("   Genotype dims for LFMM:", dim(geno_filtered), "\n")
cat("   K value:", K_lfmm, "\n")

mod_lfmm <- lfmm2(input = geno_filtered, env = X.env_scaled, K = K_lfmm, effect.sizes = TRUE)
B_mat <- mod_lfmm@B
pv <- lfmm2.test(mod_lfmm, input = geno_filtered, env = X.env_scaled, full = TRUE)
candidates <- which(-log10(pv$pvalues) > -log10(pvalue_cutoff))

cat("   Candidate SNPs identified:", length(candidates), "\n")

if (length(candidates) == 0) {
  stop("No candidate SNPs found. Try increasing pvalue_cutoff to 1e-4 or 1e-3")
}

# Create projection grid
cat("\n6. Creating projection grid...\n")
long.vec <- seq(-180, 180, length = nc)
lat.vec <- seq(-90, 90, length = nc)
coord.grid <- expand.grid(Longitude = long.vec, Latitude = lat.vec)
grid_pts <- vect(coord.grid, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
cat("   Grid dimensions:", nc, "×", nc, "=", nrow(coord.grid), "cells\n")

# Process CMIP6 models
cat("\n7. Processing CMIP6 projections...\n")
cat("   Models:", length(cmip6_models), "\n")
cat("   Scenarios:", length(ssp_scenarios), "\n")
cat("   Periods:", length(future_periods), "\n")
cat("   Total combinations:", length(cmip6_models) * length(ssp_scenarios) * length(future_periods), "\n")

# Storage for offset maps
combo_maps <- list()
all_per_model_maps <- list()

# Process each model
for (model_name in cmip6_models) {
  for (ssp in ssp_scenarios) {
    for (period in future_periods) {
      
      combo_id <- paste(ssp, period, sep = "_")
      cat("\n  Processing:", model_name, "-", combo_id, "\n")
      
      tryCatch({
        # Download future climate
        fut_raw <- cmip6_world(model = model_name, ssp = ssp, time = period,
                               var = "bioc", res = res_minutes, path = work_dir)
        fut <- to_spatraster(fut_raw)
        
        # Normalize names
        raw_names <- names(fut)
        norm_names <- normalize_bio_names(raw_names)
        names(fut) <- norm_names
        
        # Match with historical climate variables
        common_vars <- intersect(names(clim_hist), norm_names)
        if (length(common_vars) == 0) {
          cat("    No common variables, skipping\n")
          next
        }
        
        cat("    Common variables:", length(common_vars), "\n")
        
        fut_matched <- fut[[common_vars]]
        
        # Extract future environment at grid points
        Y.env_df <- terra::extract(fut_matched, grid_pts)
        Y.env <- Y.env_df[, -1, drop = FALSE]
        Y.env <- as.data.frame(lapply(Y.env, function(x) { 
          x[is.na(x)] <- mean(x, na.rm = TRUE)
          x 
        }))
        
        # Match columns with X.env
        shared_cols <- intersect(colnames(X.env), colnames(Y.env))
        if (length(shared_cols) == 0) {
          cat("    No shared columns, skipping\n")
          next
        }
        
        cat("    Shared variables:", length(shared_cols), "\n")
        
        X.env_use <- X.env[, shared_cols, drop = FALSE]
        Y.env_use <- Y.env[, shared_cols, drop = FALSE]
        
        # Scale future using historical means/sds
        means_hist <- colMeans(X.env_use)
        sds_hist <- apply(X.env_use, 2, sd)
        Y.env_scaled <- scale(Y.env_use, center = means_hist, scale = sds_hist)
        X.env_scaled_shared <- scale(X.env_use)
        
        # CRITICAL: Check B_mat dimensions and subset correctly
        # B_mat from LFMM has dimensions: SNPs × env_variables
        cat("    B_mat dimensions:", dim(B_mat), "(SNPs × env vars)\n")
        cat("    Candidate SNPs:", length(candidates), "\n")
        cat("    Shared env columns:", length(shared_cols), "\n")
        
        # Subset B matrix: use candidate SNPs (rows) and shared env vars (columns)
        B_use <- B_mat[candidates, shared_cols, drop = FALSE]
        cat("    B_use dimensions:", dim(B_use), "(candidate SNPs × shared env vars)\n")
        
        # Get genotypes for candidate SNPs only
        geno_candidates <- geno_filtered[, candidates, drop = FALSE]
        cat("    Genotype dimensions:", dim(geno_candidates), "(samples × candidate SNPs)\n")
        
        # Check dimensions before matrix multiplication
        # geno_candidates: n_samples × n_candidates
        # B_use: n_candidates × n_env_vars
        # Y.env_scaled: n_grid_cells × n_env_vars (needs transpose)
        # Result: n_samples × n_grid_cells
        
        cat("    Calculating predictions...\n")
        cat("      geno_candidates: ", dim(geno_candidates), "\n")
        cat("      B_use: ", dim(B_use), "\n")
        cat("      Y.env_scaled: ", dim(Y.env_scaled), "\n")
        
        # Genetic offset calculation
        # For future: samples × candidates %*% candidates × env_vars %*% env_vars × grid_cells
        preds_fut <- geno_candidates %*% B_use %*% t(Y.env_scaled)
        cat("      preds_fut dimensions: ", dim(preds_fut), "\n")
        
        # For current: samples × candidates %*% candidates × env_vars %*% env_vars × samples
        preds_cur <- geno_candidates %*% B_use %*% t(X.env_scaled_shared)
        cat("      preds_cur dimensions: ", dim(preds_cur), "\n")
        
        # Calculate offset: difference between future and current predictions
        # preds_fut is samples × grid_cells
        # preds_cur is samples × samples (current environment for each sample)
        # We need to expand preds_cur to match preds_fut dimensions
        
        # Take mean prediction across samples for current environment
        mean_pred_cur <- rowMeans(preds_cur)
        
        # Subtract mean current prediction from all future predictions
        offset_matrix <- abs(sweep(preds_fut, 1, mean_pred_cur, "-"))
        
        # Average offset across samples for each grid cell
        offset_per_cell <- colMeans(offset_matrix)
        
        cat("    Offset calculated for", length(offset_per_cell), "grid cells\n")
        
        # Create map matrix
        offset_map <- matrix(offset_per_cell, nrow = nc, ncol = nc, byrow = FALSE)
        
        cat("    Offset range: [", round(min(offset_map, na.rm = TRUE), 3), 
            ",", round(max(offset_map, na.rm = TRUE), 3), "]\n")
        
        # Store per-model map
        map_id <- paste(model_name, combo_id, sep = "_")
        all_per_model_maps[[map_id]] <- offset_map
        
        # Add to combo average
        if (!combo_id %in% names(combo_maps)) {
          combo_maps[[combo_id]] <- list()
        }
        combo_maps[[combo_id]][[length(combo_maps[[combo_id]]) + 1]] <- offset_map
        
        cat("    Successfully processed\n")
        
      }, error = function(e) {
        cat("    ERROR:", e$message, "\n")
      })
    }
  }
}

# Calculate combo means
cat("\n8. Calculating scenario means...\n")
combo_mean_maps <- list()
for (combo_id in names(combo_maps)) {
  maps_list <- combo_maps[[combo_id]]
  if (length(maps_list) > 0) {
    combo_mean_maps[[combo_id]] <- Reduce("+", maps_list) / length(maps_list)
    cat("   ", combo_id, ": averaged", length(maps_list), "models\n")
  }
}

# Calculate grand mean across all scenarios
cat("\n9. Calculating grand mean across all models and scenarios...\n")
valid_maps <- Filter(is_valid_map, all_per_model_maps)
cat("   Valid maps:", length(valid_maps), "\n")

if (length(valid_maps) == 0) {
  stop("No valid offset maps generated! Check climate data downloads.")
}

grand_mean <- Reduce("+", valid_maps) / length(valid_maps)
cat("   Grand mean offset range: [", round(min(grand_mean, na.rm = TRUE), 3), 
    ",", round(max(grand_mean, na.rm = TRUE), 3), "]\n")

# Determine regional extent with buffer
cat("\n10. Determining regional plotting extent...\n")
buffer_degrees <- 5
min_lon_buffered <- max(min_lon - buffer_degrees, -180)
max_lon_buffered <- min(max_lon + buffer_degrees, 180)
min_lat_buffered <- max(min_lat - buffer_degrees, -90)
max_lat_buffered <- min(max_lat + buffer_degrees, 90)

cat("   Buffered extent:\n")
cat("   Longitude: [", min_lon_buffered, ",", max_lon_buffered, "]\n")
cat("   Latitude: [", min_lat_buffered, ",", max_lat_buffered, "]\n")

# Find regional indices
lon_idx <- which(long.vec >= min_lon_buffered & long.vec <= max_lon_buffered)
lat_idx <- which(lat.vec >= min_lat_buffered & lat.vec <= max_lat_buffered)

cat("   Grid points in region: Lon =", length(lon_idx), ", Lat =", length(lat_idx), "\n")

# Extract regional data
long.vec.regional <- long.vec[lon_idx]
lat.vec.regional <- lat.vec[lat_idx]
grand_mean_regional <- grand_mean[lon_idx, lat_idx, drop = FALSE]

# Prepare color scale
global_min <- min(grand_mean_regional, na.rm = TRUE)
global_max <- max(grand_mean_regional, na.rm = TRUE)
color_fun <- colorRampPalette(c("blue", "cyan", "yellow", "orange", "red"))(100)

# Create genomic offset plot (Panel B)
cat("\n11. Creating genomic offset plot (Panel B)...\n")

# Save panel B as a PNG using base graphics
temp_panel_b <- tempfile(fileext = ".png")
png(temp_panel_b, width = 1200, height = 800, res = 150)
par(mar = c(4, 4, 3, 5))
fields::image.plot(long.vec.regional, lat.vec.regional, grand_mean_regional,
                   col = color_fun, zlim = c(global_min, global_max),
                   xlab = "Longitude", ylab = "Latitude",
                   main = "B) Genomic Offset (Grand Mean)",
                   cex.main = 1.2, font.main = 2)
points(coords$Longitude, coords$Latitude, pch = 21, bg = "white", col = "black", cex = 1.5)
maps::map(xlim = range(long.vec.regional), ylim = range(lat.vec.regional), 
          add = TRUE, interior = TRUE, col = "grey30", lwd = 0.8)
maps::map(xlim = range(long.vec.regional), ylim = range(lat.vec.regional), 
          add = TRUE, interior = FALSE, col = "grey20", lwd = 1.5)
dev.off()

cat("   Panel B created successfully\n")

# ---------------------- COMBINE PANELS ----------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("COMBINING PANELS INTO FINAL FIGURE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Convert panel B PNG to grob for combination
library(png)
library(grid)
panel_b_img <- readPNG(temp_panel_b)
panel_b_grob <- rasterGrob(panel_b_img, interpolate = TRUE)

# Create final combined plot
combined_plot <- plot_grid(
  panel_A,
  panel_b_grob,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = NULL
)

# Save combined figure
output_file <- file.path(out_dir, "Combined_RDA_GenomicOffset_Figure.png")
ggsave(output_file, combined_plot, width = 16, height = 7, dpi = 300)

# Also save individual panels
ggsave(file.path(out_dir, "Panel_A_RDA.png"), panel_A, width = 8, height = 7, dpi = 300)
file.copy(temp_panel_b, file.path(out_dir, "Panel_B_GenomicOffset.png"))

# Save summary statistics
summary_stats <- data.frame(
  Metric = c(
    "Samples analyzed",
    "SNPs (RDA)",
    "SNPs (LFMM filtered)",
    "Candidate SNPs",
    "Environmental variables (RDA)",
    "RDA variance explained (%)",
    "Climate models used",
    "Climate scenarios",
    "Min genomic offset",
    "Max genomic offset",
    "Mean genomic offset"
  ),
  Value = c(
    nrow(genotype_data),
    ncol(genotype_data),
    ncol(geno_filtered),
    length(candidates),
    ncol(env_vars_final),
    round(var_explained * 100, 2),
    length(cmip6_models),
    length(ssp_scenarios) * length(future_periods),
    round(global_min, 3),
    round(global_max, 3),
    round(mean(grand_mean_regional, na.rm = TRUE), 3)
  )
)

write.csv(summary_stats, file.path(out_dir, "Analysis_Summary.csv"), row.names = FALSE)

# Clean up temp file
unlink(temp_panel_b)

# ---------------------- FINAL SUMMARY ----------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat("Output files saved to:", out_dir, "\n\n")
cat("Main output:\n")
cat("  -", output_file, "\n\n")
cat("Summary statistics:\n")
print(summary_stats)
cat("\n")
cat("Analysis completed successfully!\n")