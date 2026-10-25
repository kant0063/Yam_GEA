################################################################################
# 07.5_EGWAS_relaxed.R                              Analysis Plan: A2 (relaxed set)
# Relaxed-set version of 07 for the missing-data comparison (script 19).
# Same model as 07 (EMMAX/P3D, Q + K for the full panel; within-species covariates),
# with changes needed for 1.57 M SNPs:
#   - SNPs scanned in chunks of 200,000 (one REML fit per trait, reused per chunk)
#   - uses the ICs saved by 07 (strict) so both sets test identical IC phenotypes
#   - writes only what 19 needs plus summaries: significant hits with effects and
#     lambda/counts (per-SNP p-values for 1.57 M SNPs x 39 traits are not saved;
#     set SAVE_ALL_P <- TRUE to keep them)
#   - no Manhattan plots
# Outputs (tables/relaxed/): TableS4_EGWAS_significant_with_effects.csv,
#   TableS4b_lambda_and_counts.csv, EGWAS_summary_by_set.csv
# Note: BH-FDR over 1.57 M tests is stricter than over 63 k; 19 therefore compares
# associated 10-kb windows, not SNP counts.
# v2: full-rank covariates (via 00), small sets skipped, NA p-values handled.
# Always resets YAM_GENO_SET to "strict" when it finishes or fails.
################################################################################
SAVE_ALL_P <- FALSE
CHUNK <- 2e5
Sys.setenv(YAM_GENO_SET = "relaxed")
tryCatch({
source("00_config_and_helpers.R"); need("rrBLUP")
stopifnot(CFG$geno_set == "relaxed")
M <- load_master(); meta <- M$meta
ic_file <- file.path("results", "strict", "ICs.rds")
if (!file.exists(ic_file)) stop("Run 07_EGWAS_revised.R (strict) first - its ICs are reused here")
IC <- readRDS(ic_file)[rownames(M$E), , drop = FALSE]
Y <- cbind(M$E, IC)

S <- subsets(M); all_hits <- list(); lam <- list()
msg("Set sizes:", paste(names(S), lengths(S), sep = " = ", collapse = "; "))
small <- names(S)[lengths(S) < CFG$min_n_within]
if (length(small)) { msg("WARNING: skipping set(s) with fewer than", CFG$min_n_within, "samples:", paste(small, collapse = ", "))
  S <- S[setdiff(names(S), small)] }
for (s in names(S)) {
  ids <- intersect(S[[s]], rownames(M$G)); m <- meta[match(ids, meta$SampleID), ]
  af_s <- colMeans(M$G[ids, ], na.rm = TRUE) / 2
  cols <- which(pmin(af_s, 1 - af_s) >= CFG$maf)
  # kinship for the subset (chunked VanRaden)
  K <- matrix(0, length(ids), length(ids)); den <- 2 * sum(af_s[cols] * (1 - af_s[cols]))
  for (st in seq(1, length(cols), by = CHUNK)) {
    j <- cols[st:min(length(cols), st + CHUNK - 1)]
    Z <- sweep(M$G[ids, j, drop = FALSE], 2, 2 * af_s[j]); Z[is.na(Z)] <- 0; K <- K + tcrossprod(Z) }
  K <- K / den; dimnames(K) <- list(ids, ids)
  X <- switch(s,
    full = M$PCs[ids, 1:CFG$n_pc],
    praehensilis = matrix(as.numeric(m$Lineage %in% "Cameroonian"), ncol = 1, dimnames = list(NULL, "Cameroonian")),
    rotundata = prcomp(M$G[ids, intersect(cols, which(M$thin))])$x[, 1:2])
  Xf <- full_rank_X(X, length(ids)); X <- Xf[, -1, drop = FALSE]; if (!ncol(X)) X <- NULL   # reduce once per set
  msg("Set", s, ":", length(ids), "samples,", length(cols), "SNPs (MAF >=", CFG$maf, ")")
  P_all <- list()
  for (tr in colnames(Y)) {
    if (!is.finite(sd(Y[ids, tr])) || sd(Y[ids, tr]) == 0) { msg("  skip", tr, "(no variation in", s, ")"); next }
    y <- if (isTRUE(CFG$int_transform)) rint(Y[ids, tr]) else as.numeric(scale(Y[ids, tr]))
    fit <- rrBLUP::mixed.solve(y, X = Xf, K = K, bounds = c(CFG$min_ve_ratio %||% 0.05, 1e9))
    res_tr <- list()
    for (st in seq(1, length(cols), by = CHUNK)) {
      j <- cols[st:min(length(cols), st + CHUNK - 1)]
      Gj <- M$G[ids, j, drop = FALSE]; storage.mode(Gj) <- "double"
      res_tr[[length(res_tr) + 1]] <- emmax_scan(y, Gj, K, X, fit = fit)
    }
    r <- bind_rows(res_tr); r$q <- p.adjust(r$p, "BH"); r$trait <- tr; r$set <- s
    lam[[paste(s, tr)]] <- data.frame(set = s, trait = tr, lambda = lambda_gc(r$p), h2_null = r$h2_null[1],
                                      n_sig = sum(r$q < CFG$fdr, na.rm = TRUE), n_tested = sum(is.finite(r$p)),
                                      boundary = r$h2_null[1] > 0.94)
    all_hits[[paste(s, tr)]] <- r %>% filter(!is.na(q), q < CFG$fdr)
    if (SAVE_ALL_P) P_all[[tr]] <- r[, c("SNP", "beta", "se", "p", "af", "PVE")]
    msg(" ", tr, ": lambda", round(lam[[paste(s, tr)]]$lambda, 2), "| hits", lam[[paste(s, tr)]]$n_sig)
  }
  if (SAVE_ALL_P) saveRDS(P_all, res(paste0("egwas_pvalues_", s, ".rds")))
}
hits <- bind_rows(all_hits)
hits <- if (nrow(hits)) left_join(hits, M$snp[, c("SNP", "CHROM", "POS")], by = "SNP") else
  data.frame(SNP = character(), beta = numeric(), se = numeric(), p = numeric(), af = numeric(), PVE = numeric(),
             q = numeric(), trait = character(), set = character(), CHROM = character(), POS = numeric())
hits$beta_original_units <- if (nrow(hits) && !isTRUE(CFG$int_transform)) hits$beta * apply(Y, 2, sd)[hits$trait] else rep(NA_real_, nrow(hits))
lam <- bind_rows(lam)
write.csv(hits, tab("TableS4_EGWAS_significant_with_effects.csv"), row.names = FALSE)
write.csv(lam, tab("TableS4b_lambda_and_counts.csv"), row.names = FALSE)
summ <- hits %>% group_by(set) %>% summarise(associations = n(), unique_markers = n_distinct(SNP),
          median_PVE_pct = 100 * median(PVE), .groups = "drop")
write.csv(summ, tab("EGWAS_summary_by_set.csv"), row.names = FALSE); print(summ)
}, finally = { Sys.setenv(YAM_GENO_SET = "strict"); message("YAM_GENO_SET reset to 'strict'") })
