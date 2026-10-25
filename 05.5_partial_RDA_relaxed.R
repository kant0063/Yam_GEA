################################################################################
# 05.5_partial_RDA_relaxed.R                        Analysis Plan: A2 (relaxed set)
# Runs the partial RDA / variance partitioning of 05 on the relaxed, imputed set.
# Output: tables/relaxed/Table2_pRDA_summary.csv (read by 19) and the other
#         Table 2 files, figures/relaxed/Fig3A_*.png.
#
# v2 - memory: the relaxed master object holds all 1.57 M SNPs (~1-2 GB) plus kinship,
# and 05 only needs the 7-kb thinned, MAF-filtered SNPs. This script
#   1. builds once a LIGHT copy (results/relaxed/master_rda.rds) containing only the
#      thinned SNPs and no kinship, and points 05 at it (YAM_MASTER);
#   2. caps the RDA response at RDA_MAX_SNPS randomly chosen thinned SNPs
#      (RDA fractions are stable well below this; set 0 for no cap);
#   3. runs permutations on RDA_CORES workers - each forked worker can duplicate
#      memory, so 1 is the safest; raise it if memory allows.
# Always resets YAM_GENO_SET / YAM_MASTER and the settings when it finishes or fails.
################################################################################
RDA_MAX_SNPS <- 50000
RDA_CORES    <- 1
Sys.setenv(YAM_GENO_SET = "relaxed")
tryCatch({
  full_master <- file.path("results", "relaxed", "master.rds")
  light       <- file.path("results", "relaxed", "master_rda.rds")
  if (!file.exists(full_master)) stop("Run 01.5_build_master_relaxed.R first")
  if (!file.exists(light) || file.mtime(light) < file.mtime(full_master)) {
    message("Building light master for RDA ...")
    M <- readRDS(full_master)
    af <- colMeans(M$G, na.rm = TRUE) / 2
    keep <- which(M$thin & pmin(af, 1 - af) >= 0.01)         # thinned + not monomorphic (05 re-filters MAF per set)
    M$G <- M$G[, keep, drop = FALSE]; M$snp <- M$snp[keep, ]; M$thin <- rep(TRUE, length(keep)); M$K <- NULL
    saveRDS(M, light); rm(M); invisible(gc())
    message("Light master: ", length(keep), " thinned SNPs -> ", light)
  }
  Sys.setenv(YAM_MASTER = light, YAM_RDA_MAX_SNPS = RDA_MAX_SNPS, YAM_RDA_CORES = RDA_CORES)
  source("05_partial_RDA_varpart.R", echo = FALSE)
  msg("Relaxed pRDA done:", file.path("tables", "relaxed", "Table2_pRDA_summary.csv"))
}, finally = {
  Sys.setenv(YAM_GENO_SET = "strict"); Sys.unsetenv(c("YAM_MASTER", "YAM_RDA_MAX_SNPS", "YAM_RDA_CORES"))
  message("YAM_GENO_SET reset to 'strict'; RDA settings cleared")
})
