################################################################################
# 01.5_build_master_relaxed.R                       Analysis Plan: A2 (relaxed set)
# Relaxed-set version of 01 for the missing-data sensitivity analysis (script 19).
# Run AFTER 01 (strict) and 03 (creates Yam_relaxed_imputed.vcf.gz).
#
# Differences from 01, needed for 1.57 M SNPs:
#   - reads genotypes with bcftools query -> integer dosage (no vcfR; far less memory)
#   - takes metadata, environment, VIF variables and sample order from the STRICT
#     master object, so both sets use identical samples and predictors
#   - kinship computed in chunks (VanRaden), PCs from MAF-filtered SNPs thinned to 7 kb
# Output: results/relaxed/master.rds (same structure as 01, read by 05.5, 07.5, 19)
# Memory: ~3-6 GB. Always resets YAM_GENO_SET to "strict" when it finishes or fails.
################################################################################
Sys.setenv(YAM_GENO_SET = "relaxed")
tryCatch({
source("00_config_and_helpers.R")
stopifnot(CFG$geno_set == "relaxed")
strict_master <- file.path("results", "strict", "master.rds")
if (!file.exists(strict_master)) stop("Run 01_build_master_data.R (strict) first")
if (!file.exists(CFG$vcf_relaxed)) stop("Run 03_relaxed_missing_dataset.R first: ", CFG$vcf_relaxed, " not found")
bcftools <- Sys.getenv("YAM_BCFTOOLS", Sys.which("bcftools")); if (!nzchar(bcftools)) bcftools <- "/opt/homebrew/bin/bcftools"
MS <- readRDS(strict_master); ids <- MS$meta$SampleID

# ---- genotypes: bcftools query -> 0/1/2 integers ---------------------------------
vcf_samples <- system(paste(shQuote(bcftools), "query -l", shQuote(CFG$vcf_relaxed)), intern = TRUE)
keep_s <- intersect(ids, vcf_samples)
if (length(keep_s) < length(ids)) msg("WARNING:", length(ids) - length(keep_s), "strict samples missing from relaxed VCF")
writeLines(keep_s, res("samples_keep.txt"))
msg("Reading", CFG$vcf_relaxed, "for", length(keep_s), "samples ...")
cmd <- sprintf("%s query -S %s -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' %s | sed -e 's/[|]/\\//g' -e 's/0\\/0/0/g' -e 's/0\\/1/1/g' -e 's/1\\/0/1/g' -e 's/1\\/1/2/g' -e 's/\\.\\/\\./NA/g'",
               shQuote(bcftools), shQuote(res("samples_keep.txt")), shQuote(CFG$vcf_relaxed))
dt <- fread(cmd = cmd, header = FALSE, sep = "\t", na.strings = "NA")
snp <- data.frame(SNP = paste(dt$V1, dt$V2, sep = "_"), CHROM = dt$V1, POS = as.numeric(dt$V2),
                  REF = dt$V3, ALT = dt$V4, stringsAsFactors = FALSE)
G <- base::t(as.matrix(dt[, -(1:4)])); rm(dt); gc()
storage.mode(G) <- "integer"; rownames(G) <- keep_s; colnames(G) <- snp$SNP
msg("Relaxed genotypes:", nrow(G), "samples x", ncol(G), "SNPs; missing =", sum(is.na(G)))
rg <- range(G, na.rm = TRUE); if (rg[1] < 0 || rg[2] > 2) stop("Unexpected genotype codes - check the sed conversion")

# ---- reuse strict metadata / environment / predictors ------------------------------
meta <- MS$meta[match(keep_s, MS$meta$SampleID), ]
E <- MS$E[keep_s, , drop = FALSE]

# ---- PCs from MAF-filtered, 7-kb thinned SNPs ---------------------------------------
af <- colMeans(G, na.rm = TRUE) / 2; maf_ok <- pmin(af, 1 - af) >= CFG$maf
thin <- thin_snps(snp)
sel <- which(maf_ok & thin)
pca <- prcomp(G[, sel], center = TRUE, scale. = FALSE)
PCs <- pca$x[, 1:10]; pc_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
msg("PCA on", length(sel), "thinned SNPs; PC1-3 %:", paste(round(pc_var[1:3], 2), collapse = ", "))

# ---- kinship (VanRaden 2008), chunked to limit memory --------------------------------
p <- af; den <- 2 * sum(p * (1 - p)); K <- matrix(0, nrow(G), nrow(G), dimnames = list(keep_s, keep_s))
for (st in seq(1, ncol(G), by = 2e5)) {
  j <- st:min(ncol(G), st + 2e5 - 1)
  Z <- sweep(G[, j, drop = FALSE], 2, 2 * p[j]); Z[is.na(Z)] <- 0
  K <- K + tcrossprod(Z)
}
K <- K / den

M <- list(G = G, snp = snp, meta = meta, E = E, env_vars = colnames(E),
          vif_keep = MS$vif_keep, vif_keep_bio = MS$vif_keep_bio, PCs = PCs, pc_var = pc_var,
          K = K, thin = thin, cfg = CFG)
saveRDS(M, res("master.rds"))
msg("Saved", res("master.rds"), "| SNPs:", ncol(G), "| thinned for RDA/PCA:", sum(thin))
}, finally = { Sys.setenv(YAM_GENO_SET = "strict"); message("YAM_GENO_SET reset to 'strict'") })
