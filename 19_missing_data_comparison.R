################################################################################
# 19_missing_data_comparison.R                             Analysis Plan: A2 (part 2)
# Compares strict (zero-missing) and relaxed (imputed) results after running
# 01, 05 and 07 for both sets. Output: Table S7b.               New script.
################################################################################
Sys.setenv(YAM_GENO_SET = "strict"); source("00_config_and_helpers.R"); need("vegan")
rpath <- function(set, f) file.path("results", set, f); tpath <- function(set, f) file.path("tables", set, f)   # (named so they never mask base::t)
Ms <- readRDS(rpath("strict", "master.rds")); Mr <- readRDS(rpath("relaxed", "master.rds"))
ids <- intersect(rownames(Ms$PCs), rownames(Mr$PCs))
pro <- vegan::protest(Ms$PCs[ids, 1:4], Mr$PCs[ids, 1:4], permutations = 999)
rs <- read.csv(tpath("strict", "Table2_pRDA_summary.csv")); rr <- read.csv(tpath("relaxed", "Table2_pRDA_summary.csv"))
hs <- read.csv(tpath("strict", "TableS4_EGWAS_significant_with_effects.csv")); hr <- read.csv(tpath("relaxed", "TableS4_EGWAS_significant_with_effects.csv"))
win <- function(h) paste(h$set, h$trait, h$CHROM, floor(h$POS / CFG$wza_window))     # compare at window level (SNP sets differ)
s7b <- data.frame(
  metric = c("SNPs", "Procrustes r (PC1-4)", "Procrustes P", "pRDA env|anc+geo adj R2 (full)",
             "EGWAS associations (full)", "EGWAS windows shared (all sets)"),
  strict = c(ncol(Ms$G), NA, NA, rs$env_given_geo_anc_adjR2[rs$set == "full"], sum(hs$set == "full"), NA),
  relaxed = c(ncol(Mr$G), pro$t0, pro$signif, rr$env_given_geo_anc_adjR2[rr$set == "full"], sum(hr$set == "full"),
              length(intersect(win(hs), win(hr)))))
dir.create("tables/comparison", showWarnings = FALSE, recursive = TRUE)
write.csv(s7b, "tables/comparison/TableS7b_strict_vs_relaxed.csv", row.names = FALSE); print(s7b)
