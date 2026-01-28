# =============================================================================
# INTEGRATED POPULATION GENOMICS ANALYSIS
# Combines: ROH Analysis, Private Alleles, and FST Outlier Detection
# =============================================================================
# This script performs comprehensive population genomics analysis including:
# 1. Runs of Homozygosity (ROH) detection and Froh calculation
# 2. Private allele identification by population
# 3. FST outlier analysis
# 4. Creates integrated three-panel publication figure
#
# Author: Created for Michael Kantar
# Date: 2026-01-27
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(detectRUNS)
  library(data.table)
  library(tidyverse)
  library(vcfR)
  library(hierfstat)
  library(OutFLANK)
  library(qvalue)
  library(MetBrewer)
  library(patchwork)
})

# =============================================================================
# SECTION 1: PARAMETERS - CUSTOMIZE THESE
# =============================================================================

cat("=============================================================================\n")
cat("INTEGRATED POPULATION GENOMICS ANALYSIS\n")
cat("=============================================================================\n\n")

# --- File paths ---
# PLINK files for ROH analysis
ped_file <- "Yam_no_missing_plink.ped"
map_file <- "Yam_no_missing_plink.map"

# VCF files
vcf_private_alleles <- "Calling_ALL_Rotundata_Allc05.vcf.gz"  # For private alleles
vcf_fst <- "Yam_no_missing.vcf"  # For FST analysis

# Population metadata files
pop_file_roh <- "Yam_runs_of_homozygozity_meta.csv"  # For ROH
pop_file_private <- "Yam_private_allele_script.csv"  # For private alleles
pop_file_private_detailed <- "Yam_private_allele_script_pop.csv"  # Detailed pop for private alleles
pop_file_fst <- "Yam_meta_data.csv"  # For FST

# --- ROH parameters ---
min_snps_roh <- 20
max_het <- 1
max_opp <- 1
max_miss <- 1

# --- FST outlier thresholds ---
fst_threshold_995 <- 0.995
fst_threshold_99 <- 0.99
fst_threshold_95 <- 0.95

# =============================================================================
# SECTION 2: ROH ANALYSIS
# =============================================================================

cat("\n=============================================================================\n")
cat("PART 1: RUNS OF HOMOZYGOSITY (ROH) ANALYSIS\n")
cat("=============================================================================\n\n")

# --- 2.1: Standardize PLINK files ---
cat("Standardizing PLINK files for detectRUNS...\n")

# Read and standardize MAP file
map_data <- fread(map_file, header = FALSE, sep = "\t")
if (ncol(map_data) == 1) {
  map_data <- fread(map_file, header = FALSE, sep = " ")
}
if (ncol(map_data) != 4) {
  stop(sprintf("MAP file has %d columns, expected 4", ncol(map_data)))
}
colnames(map_data) <- c("CHR", "SNP", "cM", "POS")
cat(sprintf("  ✓ MAP file: %d SNPs\n", nrow(map_data)))

# Write standardized MAP
write.table(map_data, "temp_standard.map", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Read and standardize PED file
ped_data <- fread(ped_file, header = FALSE)
write.table(ped_data, "temp_standard.ped", 
            quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
cat(sprintf("  ✓ PED file standardized\n\n"))

rm(ped_data)
gc()

# --- 2.2: Run ROH detection ---
cat("Running ROH detection (consecutive method)...\n")

roh_consecutive <- consecutiveRUNS.run(
  genotypeFile = "temp_standard.ped",
  mapFile = "temp_standard.map",
  minSNP = min_snps_roh,
  ROHet = TRUE,
  maxOppRun = max_opp,
  maxMissRun = max_miss
)

cat("  ✓ ROH detection completed\n\n")

# --- 2.3: Merge with population data ---
pop_data_roh <- read.csv(pop_file_roh, stringsAsFactors = FALSE)
roh_results <- roh_consecutive %>%
  left_join(pop_data_roh, by = c("id" = "sample_id"))

cat(sprintf("  ✓ Total ROH detected: %d\n", nrow(roh_results)))

# Check the structure of ROH results for debugging
cat("\nROH results structure:\n")
cat("Columns:", paste(colnames(roh_results), collapse = ", "), "\n")
cat("First few rows:\n")
print(head(roh_results, 3))
cat("\n")

# --- 2.4: Calculate Froh statistics ---
cat("Calculating Froh coefficients...\n")

# Calculate genome length
genome_length <- map_data %>% 
  group_by(CHR) %>% 
  summarise(max_pos = max(POS)) %>% 
  pull(max_pos) %>%
  sum()

# Calculate Froh for each individual
# detectRUNS returns 'lengthBps' column
froh_values <- roh_results %>%
  filter(!is.na(population)) %>%
  group_by(id, population) %>%
  summarise(
    total_roh_length = sum(lengthBps),
    n_roh = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Froh = total_roh_length / genome_length
  )

# Add individuals with no ROH
all_samples <- pop_data_roh %>%
  filter(!is.na(population))

samples_with_roh <- unique(froh_values$id)
samples_without_roh <- setdiff(all_samples$sample_id, samples_with_roh)

if (length(samples_without_roh) > 0) {
  zero_roh <- data.frame(
    id = samples_without_roh,
    population = all_samples$population[match(samples_without_roh, all_samples$sample_id)],
    total_roh_length = 0,
    n_roh = 0,
    Froh = 0
  )
  froh_values <- bind_rows(froh_values, zero_roh)
}

cat(sprintf("  ✓ Froh calculated for %d individuals\n\n", nrow(froh_values)))

# Summary by population
froh_summary_pop <- froh_values %>%
  group_by(population) %>%
  summarise(
    n_individuals = n(),
    mean_Froh = mean(Froh),
    median_Froh = median(Froh),
    sd_Froh = sd(Froh),
    min_Froh = min(Froh),
    max_Froh = max(Froh),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_Froh))

cat("Froh summary by population:\n")
print(froh_summary_pop)
cat("\n")

# =============================================================================
# SECTION 3: PRIVATE ALLELE ANALYSIS
# =============================================================================

cat("=============================================================================\n")
cat("PART 2: PRIVATE ALLELE ANALYSIS\n")
cat("=============================================================================\n\n")

# --- 3.1: Load VCF and population data ---
cat("Reading VCF file for private alleles...\n")
vcf_pa <- read.vcfR(vcf_private_alleles)
cat(sprintf("  ✓ VCF loaded: %d variants, %d samples\n", nrow(vcf_pa@fix), ncol(vcf_pa@gt) - 1))

pop_data_pa <- read.csv(pop_file_private, stringsAsFactors = FALSE)
cat(sprintf("  ✓ Population data loaded: %d samples\n\n", nrow(pop_data_pa)))

# --- 3.2: Extract genotypes ---
gt_pa <- extract.gt(vcf_pa, element = "GT")
vcf_samples <- colnames(gt_pa)

# Ensure population data matches VCF
pop_data_pa <- pop_data_pa %>%
  filter(sample_id %in% vcf_samples) %>%
  arrange(match(sample_id, vcf_samples))

# --- 3.3: Count private alleles ---
cat("Counting private alleles by population...\n")

extract_alleles <- function(genotype) {
  if (is.na(genotype) || genotype == "./.") {
    return(character(0))
  }
  alleles <- unlist(strsplit(genotype, "[/|]"))
  alleles <- alleles[alleles != "."]
  return(unique(alleles))
}

populations <- unique(pop_data_pa$population)
private_allele_counts <- data.frame(
  population = populations,
  private_alleles = 0,
  stringsAsFactors = FALSE
)

# Progress bar
pb <- txtProgressBar(min = 0, max = nrow(gt_pa), style = 3)

for (i in 1:nrow(gt_pa)) {
  setTxtProgressBar(pb, i)
  
  # Get alleles in each population
  pop_alleles <- list()
  for (pop in populations) {
    pop_samples <- pop_data_pa$sample_id[pop_data_pa$population == pop]
    pop_samples <- pop_samples[pop_samples %in% vcf_samples]
    pop_gts <- gt_pa[i, pop_samples]
    alleles <- unique(unlist(lapply(pop_gts, extract_alleles)))
    pop_alleles[[pop]] <- alleles
  }
  
  # Find private alleles
  for (pop in populations) {
    if (length(pop_alleles[[pop]]) > 0) {
      other_pops <- setdiff(populations, pop)
      other_alleles <- unique(unlist(pop_alleles[other_pops]))
      private <- setdiff(pop_alleles[[pop]], other_alleles)
      
      if (length(private) > 0) {
        private_allele_counts$private_alleles[private_allele_counts$population == pop] <- 
          private_allele_counts$private_alleles[private_allele_counts$population == pop] + length(private)
      }
    }
  }
}
close(pb)

cat("\n  ✓ Private allele counting completed\n\n")

# Summary statistics
private_summary <- pop_data_pa %>%
  group_by(population) %>%
  summarise(n_samples = n(), .groups = 'drop') %>%
  left_join(private_allele_counts, by = "population") %>%
  mutate(
    private_alleles_per_sample = round(private_alleles / n_samples, 2),
    percent_of_total = round(100 * private_alleles / sum(private_alleles), 2)
  ) %>%
  arrange(desc(private_alleles))

cat("Private allele summary:\n")
print(private_summary)
cat("\n")

# =============================================================================
# SECTION 4: FST OUTLIER ANALYSIS
# =============================================================================

cat("=============================================================================\n")
cat("PART 3: FST OUTLIER ANALYSIS\n")
cat("=============================================================================\n\n")

# --- 4.1: Load VCF for FST ---
cat("Reading VCF file for FST analysis...\n")
vcf_fst_data <- read.vcfR(vcf_fst)
cat(sprintf("  ✓ VCF loaded: %d SNPs, %d samples\n", nrow(vcf_fst_data@fix), ncol(vcf_fst_data@gt) - 1))

# Extract genotypes
gt_numeric_fst <- extract.gt(vcf_fst_data, element = "GT", as.numeric = TRUE)
genotype_matrix_fst <- t(gt_numeric_fst)

# Get SNP information
snp_info_fst <- data.frame(
  SNP_index = 1:nrow(vcf_fst_data@fix),
  CHROM = vcf_fst_data@fix[, "CHROM"],
  POS = as.numeric(vcf_fst_data@fix[, "POS"]),
  stringsAsFactors = FALSE
)

# --- 4.2: Load population metadata ---
pop_data_fst <- read.csv(pop_file_fst, stringsAsFactors = FALSE)
cat(sprintf("  ✓ Metadata loaded: %d samples\n", nrow(pop_data_fst)))

# --- 4.3: Match samples and prepare data ---
common_samples_fst <- intersect(rownames(genotype_matrix_fst), pop_data_fst$SampleID)
genotype_matched_fst <- genotype_matrix_fst[common_samples_fst, ]
pop_matched_fst <- pop_data_fst[pop_data_fst$SampleID %in% common_samples_fst, ]
pop_matched_fst <- pop_matched_fst[match(common_samples_fst, pop_matched_fst$SampleID), ]

cat(sprintf("  ✓ Matched %d samples\n", length(common_samples_fst)))

# Impute missing data
genotype_final_fst <- apply(genotype_matched_fst, 2, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
})

# --- 4.4: Calculate FST using hierfstat ---
cat("\nCalculating population FST values...\n")

# Prepare data for hierfstat
hierfstat_data <- data.frame(
  population = as.integer(factor(pop_matched_fst$Population.assignment)),
  genotype_final_fst
)

# Calculate FST per locus
fst_per_locus <- apply(hierfstat_data[, -1], 2, function(snp) {
  temp_data <- data.frame(pop = hierfstat_data$population, snp = snp)
  temp_data <- temp_data[!is.na(temp_data$snp), ]
  if (length(unique(temp_data$pop)) < 2 || length(unique(temp_data$snp)) < 2) {
    return(NA)
  }
  tryCatch({
    fst_val <- betas(temp_data, nboot = 0)$betaiovl
    return(fst_val)
  }, error = function(e) {
    return(NA)
  })
})

# Create FST results
fst_results <- snp_info_fst %>%
  mutate(FST = fst_per_locus)

cat(sprintf("  ✓ FST calculated for %d SNPs\n", sum(!is.na(fst_results$FST))))

# --- 4.5: Identify outliers ---
cat("\nIdentifying FST outliers...\n")

# Calculate quantile thresholds
threshold_995 <- quantile(fst_results$FST, probs = fst_threshold_995, na.rm = TRUE)
threshold_99 <- quantile(fst_results$FST, probs = fst_threshold_99, na.rm = TRUE)
threshold_95 <- quantile(fst_results$FST, probs = fst_threshold_95, na.rm = TRUE)

fst_results <- fst_results %>%
  mutate(
    Outlier_995 = FST >= threshold_995,
    Outlier_99 = FST >= threshold_99,
    Outlier_95 = FST >= threshold_95
  )

cat(sprintf("  ✓ FST outliers (99.5%%): %d (%.2f%%)\n", 
            sum(fst_results$Outlier_995, na.rm = TRUE),
            100 * sum(fst_results$Outlier_995, na.rm = TRUE) / nrow(fst_results)))
cat(sprintf("  ✓ FST outliers (99%%): %d (%.2f%%)\n", 
            sum(fst_results$Outlier_99, na.rm = TRUE),
            100 * sum(fst_results$Outlier_99, na.rm = TRUE) / nrow(fst_results)))
cat(sprintf("  ✓ FST outliers (95%%): %d (%.2f%%)\n\n", 
            sum(fst_results$Outlier_95, na.rm = TRUE),
            100 * sum(fst_results$Outlier_95, na.rm = TRUE) / nrow(fst_results)))

# --- 4.6: Identify private alleles in FST dataset ---
cat("Identifying private alleles in FST dataset...\n")

# Load detailed population file
pop_data_pa_detailed <- read.csv(pop_file_private_detailed, stringsAsFactors = FALSE)
colnames(pop_data_pa_detailed) <- c("SampleID", "Population")
pop_data_pa_detailed$SampleID <- gsub("\\.bam$", "", pop_data_pa_detailed$SampleID)

# Get pure populations (exclude hybrids)
all_pops <- unique(pop_data_pa_detailed$Population)
pure_pops <- all_pops[!grepl(" x ", all_pops)]

cat(sprintf("  Pure populations: %s\n", paste(pure_pops, collapse = ", ")))

# Match samples
vcf_samples_fst <- colnames(vcf_fst_data@gt)[-1]
vcf_samples_clean <- gsub("\\.bam$", "", vcf_samples_fst)
common_pa <- intersect(vcf_samples_clean, pop_data_pa_detailed$SampleID)
pop_assignments_pa <- pop_data_pa_detailed$Population[match(common_pa, pop_data_pa_detailed$SampleID)]

cat(sprintf("  Matched %d samples for private allele detection\n", length(common_pa)))

# Initialize private allele tracking
private_allele_snps <- data.frame(
  SNP_index = 1:nrow(vcf_fst_data@fix),
  stringsAsFactors = FALSE
)

for (pop in pure_pops) {
  private_allele_snps[[paste0("Private_", gsub(" ", "_", pop))]] <- FALSE
}

# Extract genotypes for private allele detection
gt_matrix_pa <- extract.gt(vcf_fst_data, element = "GT")

# Progress bar
cat("  Processing SNPs for private allele detection...\n")
pb2 <- txtProgressBar(min = 0, max = nrow(gt_matrix_pa), style = 3)

for (i in 1:nrow(gt_matrix_pa)) {
  setTxtProgressBar(pb2, i)
  
  snp_gts <- gt_matrix_pa[i, ]
  
  # Convert to numeric
  snp_numeric <- sapply(snp_gts, function(gt) {
    if (is.na(gt) || gt == "./.") return(NA)
    alleles <- strsplit(gt, "[/|]")[[1]]
    sum(as.numeric(alleles))
  })
  
  names(snp_numeric) <- gsub("\\.bam$", "", names(snp_numeric))
  snp_matched <- snp_numeric[common_pa]
  
  # Check each population for private alleles
  for (pop in pure_pops) {
    pop_samples <- common_pa[pop_assignments_pa == pop]
    other_samples <- common_pa[pop_assignments_pa != pop & !grepl(" x ", pop_assignments_pa)]
    
    if (length(pop_samples) > 0 && length(other_samples) > 0) {
      pop_gts <- snp_matched[pop_samples]
      other_gts <- snp_matched[other_samples]
      
      # Check for private alleles
      pop_has_alt <- any(pop_gts > 0, na.rm = TRUE)
      others_have_alt <- any(other_gts > 0, na.rm = TRUE)
      pop_has_ref <- any(pop_gts < 2, na.rm = TRUE)
      others_have_ref <- any(other_gts < 2, na.rm = TRUE)
      
      is_private <- (pop_has_alt && !others_have_alt) || (pop_has_ref && !others_have_ref)
      private_allele_snps[i, paste0("Private_", gsub(" ", "_", pop))] <- is_private
    }
  }
}
close(pb2)

# Add summary columns
private_cols_snp <- grep("^Private_", colnames(private_allele_snps), value = TRUE)
private_allele_snps$Is_Private_Any <- rowSums(private_allele_snps[, private_cols_snp]) > 0
private_allele_snps$Private_To_Population <- apply(private_allele_snps[, private_cols_snp], 1, function(x) {
  if (sum(x) == 0) return("None")
  pops <- gsub("Private_", "", names(x)[x])
  pops <- gsub("_", " ", pops)
  paste(pops, collapse = "; ")
})

cat("\n  ✓ Private allele detection completed\n")
cat(sprintf("  ✓ SNPs with private alleles: %d (%.2f%%)\n\n",
            sum(private_allele_snps$Is_Private_Any),
            100 * sum(private_allele_snps$Is_Private_Any) / nrow(private_allele_snps)))

# --- 4.7: Combine FST and private allele results ---
combined_results <- fst_results %>%
  left_join(private_allele_snps, by = "SNP_index")

# =============================================================================
# SECTION 5: CREATE INTEGRATED THREE-PANEL FIGURE
# =============================================================================

cat("=============================================================================\n")
cat("CREATING INTEGRATED THREE-PANEL FIGURE\n")
cat("=============================================================================\n\n")

# --- PANEL A: Private Alleles by Population ---
cat("Creating Panel A: Private alleles by population...\n")

panel_a <- ggplot(private_summary, 
                  aes(x = reorder(population, -private_alleles), 
                      y = private_alleles, 
                      fill = population)) +
  geom_col(color = "black", size = 0.3) +
  geom_text(aes(label = private_alleles), 
            vjust = -0.5, size = 3) +
  labs(
    title = "A) Private Alleles by Population",
    x = "Population",
    y = "Number of Private Alleles"
  ) +
  scale_fill_viridis_d(option = "viridis") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )

# --- PANEL B: Froh by Population ---
cat("Creating Panel B: Froh by population...\n")

panel_b <- ggplot(froh_values, 
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
  scale_fill_viridis_d(option = "plasma") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )

# --- PANEL C: Manhattan Plot of FST with Private Alleles ---
cat("Creating Panel C: Manhattan plot of FST with private alleles...\n")

# Calculate cumulative positions
chr_lengths <- snp_info_fst %>%
  group_by(CHROM) %>%
  summarise(chr_length = max(POS), .groups = "drop")

chr_lengths_cumul <- chr_lengths %>%
  arrange(CHROM) %>%
  mutate(
    chr_start = cumsum(lag(chr_length, default = 0)),
    chr_mid = chr_start + chr_length / 2,
    chr_end = chr_start + chr_length
  )

# Add cumulative positions to results
manhattan_data <- combined_results %>%
  left_join(chr_lengths_cumul %>% select(CHROM, chr_start), by = "CHROM") %>%
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

# Create Manhattan plot
panel_c <- ggplot(manhattan_data, aes(x = pos_cumul / 1e6, y = FST)) +
  # Background: alternating chromosome colors
  geom_rect(data = chr_lengths_cumul %>% mutate(odd = row_number() %% 2 == 1),
            aes(xmin = chr_start / 1e6, xmax = chr_end / 1e6,
                ymin = -Inf, ymax = Inf, fill = odd),
            inherit.aes = FALSE, alpha = 0.15) +
  scale_fill_manual(values = c("TRUE" = "gray95", "FALSE" = "white"), guide = "none") +
  # All neutral SNPs (gray)
  geom_point(data = manhattan_data %>% filter(point_type == "Neutral"),
             aes(color = point_type), size = 0.5, alpha = 0.3) +
  # Private alleles only (green)
  geom_point(data = manhattan_data %>% filter(point_type == "Private allele only"),
             aes(color = point_type), size = 1.2, alpha = 0.8) +
  # FST outliers only (orange)
  geom_point(data = manhattan_data %>% filter(point_type == "FST outlier only"),
             aes(color = point_type), size = 1.2, alpha = 0.8) +
  # Both FST outlier and private (red)
  geom_point(data = manhattan_data %>% filter(point_type == "FST outlier + Private"),
             aes(color = point_type), size = 1.5, alpha = 1) +
  # Threshold lines
  geom_hline(yintercept = threshold_995, linetype = "dashed", 
             color = "red", size = 0.5) +
  geom_hline(yintercept = threshold_99, linetype = "dashed", 
             color = "orange", size = 0.4) +
  # Color scale
  scale_color_manual(
    values = c(
      "Neutral" = "gray60",
      "Private allele only" = "#2CA02C",
      "FST outlier only" = "#FF7F0E",
      "FST outlier + Private" = "#D62728"
    ),
    name = "SNP Type"
  ) +
  # Chromosome labels
  scale_x_continuous(
    breaks = chr_lengths_cumul$chr_mid / 1e6,
    labels = chr_lengths_cumul$CHROM,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "C) FST Manhattan Plot with Private Alleles",
    x = "Chromosome",
    y = expression(F[ST])
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.border = element_rect(fill = NA, color = "gray50")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# --- Combine panels ---
cat("Combining panels into final figure...\n")

combined_figure <- (panel_a | panel_b) / panel_c +
  plot_layout(heights = c(1, 1.2))

# Save figure
ggsave("integrated_population_genomics_figure.png", 
       combined_figure, 
       width = 14, height = 10, dpi = 300)

cat("  ✓ Figure saved: integrated_population_genomics_figure.png\n\n")

# Save high-resolution version
ggsave("integrated_population_genomics_figure_hires.png", 
       combined_figure, 
       width = 14, height = 10, dpi = 600)

cat("  ✓ High-res figure saved: integrated_population_genomics_figure_hires.png\n\n")

# =============================================================================
# SECTION 6: SAVE RESULTS
# =============================================================================

cat("=============================================================================\n")
cat("SAVING ANALYSIS RESULTS\n")
cat("=============================================================================\n\n")

# ROH results
write.csv(roh_results, "roh_all_runs.csv", row.names = FALSE)
write.csv(froh_values, "froh_coefficients.csv", row.names = FALSE)
write.csv(froh_summary_pop, "froh_summary_population.csv", row.names = FALSE)

# Private allele results
write.csv(private_allele_counts, "private_allele_counts.csv", row.names = FALSE)
write.csv(private_summary, "private_allele_summary.csv", row.names = FALSE)

# FST results
write.csv(fst_results, "fst_results.csv", row.names = FALSE)

# Combined results
write.csv(combined_results, "integrated_fst_private_allele_results.csv", row.names = FALSE)

# FST outliers that are also private alleles
fst_private_outliers <- combined_results %>%
  filter(Outlier_995 & Is_Private_Any) %>%
  arrange(desc(FST))

write.csv(fst_private_outliers, "fst_outliers_with_private_alleles.csv", row.names = FALSE)

cat("CSV files saved:\n")
cat("  - roh_all_runs.csv\n")
cat("  - froh_coefficients.csv\n")
cat("  - froh_summary_population.csv\n")
cat("  - private_allele_counts.csv\n")
cat("  - private_allele_summary.csv\n")
cat("  - fst_results.csv\n")
cat("  - integrated_fst_private_allele_results.csv\n")
cat("  - fst_outliers_with_private_alleles.csv\n\n")

# =============================================================================
# SECTION 7: SUMMARY STATISTICS
# =============================================================================

cat("=============================================================================\n")
cat("INTEGRATED ANALYSIS SUMMARY\n")
cat("=============================================================================\n\n")

cat("ROH ANALYSIS:\n")
cat(sprintf("  Total ROH detected: %d\n", nrow(roh_results)))
cat(sprintf("  Individuals analyzed: %d\n", nrow(froh_values)))
cat(sprintf("  Mean Froh across all individuals: %.4f\n", mean(froh_values$Froh)))
cat(sprintf("  Populations: %d\n", length(unique(froh_values$population))))
cat("\n")

cat("PRIVATE ALLELE ANALYSIS:\n")
cat(sprintf("  Total private alleles detected: %d\n", sum(private_allele_counts$private_alleles)))
cat(sprintf("  Mean private alleles per population: %.1f\n", mean(private_summary$private_alleles)))
cat(sprintf("  Range: %d - %d\n", min(private_summary$private_alleles), max(private_summary$private_alleles)))
cat(sprintf("  Populations: %d\n", nrow(private_summary)))
cat("\n")

cat("FST OUTLIER ANALYSIS:\n")
cat(sprintf("  Total SNPs analyzed: %d\n", nrow(fst_results)))
cat(sprintf("  FST outliers (99.5%%): %d (%.2f%%)\n", 
            sum(fst_results$Outlier_995, na.rm = TRUE),
            100 * sum(fst_results$Outlier_995, na.rm = TRUE) / nrow(fst_results)))
cat(sprintf("  SNPs with private alleles: %d (%.2f%%)\n",
            sum(combined_results$Is_Private_Any, na.rm = TRUE),
            100 * sum(combined_results$Is_Private_Any, na.rm = TRUE) / nrow(combined_results)))
cat(sprintf("  SNPs that are BOTH FST outliers AND private: %d (%.2f%%)\n",
            sum(combined_results$Outlier_995 & combined_results$Is_Private_Any, na.rm = TRUE),
            100 * sum(combined_results$Outlier_995 & combined_results$Is_Private_Any, na.rm = TRUE) / nrow(combined_results)))
cat("\n")

cat("TOP 5 POPULATIONS BY FROH:\n")
top_froh <- head(froh_summary_pop, 5)
print(top_froh[, c("population", "n_individuals", "mean_Froh")], row.names = FALSE)
cat("\n")

cat("TOP 5 POPULATIONS BY PRIVATE ALLELES:\n")
top_private <- head(private_summary, 5)
print(top_private[, c("population", "private_alleles", "private_alleles_per_sample")], row.names = FALSE)
cat("\n")

if (nrow(fst_private_outliers) > 0) {
  cat("TOP 10 SNPs (FST OUTLIERS + PRIVATE ALLELES):\n")
  top_combined <- head(fst_private_outliers %>% 
                         select(SNP_index, CHROM, POS, FST, Private_To_Population), 10)
  print(top_combined, row.names = FALSE)
  cat("\n")
}

# =============================================================================
# SECTION 8: CLEANUP
# =============================================================================

cat("=============================================================================\n")
cat("CLEANING UP\n")
cat("=============================================================================\n\n")

if (file.exists("temp_standard.ped")) file.remove("temp_standard.ped")
if (file.exists("temp_standard.map")) file.remove("temp_standard.map")

cat("  ✓ Temporary files removed\n\n")

cat("=============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=============================================================================\n\n")

cat("Main output file:\n")
cat("  ★ integrated_population_genomics_figure.png\n")
cat("  ★ integrated_population_genomics_figure_hires.png (600 DPI)\n\n")

cat("This figure shows:\n")
cat("  Panel A: Number of private alleles in each population\n")
cat("  Panel B: Distribution of inbreeding coefficients (Froh) by population\n")
cat("  Panel C: Manhattan plot of FST values with private alleles highlighted\n\n")

cat("Interpretation:\n")
cat("  - Red points: SNPs under selection (high FST) that are also private to populations\n")
cat("  - Orange points: SNPs under selection but not private\n")
cat("  - Green points: Private alleles not under strong selection\n")
cat("  - Gray points: Neutral variation\n\n")

cat("=============================================================================\n")