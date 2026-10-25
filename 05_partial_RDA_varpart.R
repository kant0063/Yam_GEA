################################################################################
# 05_partial_RDA_varpart.R                                 Analysis Plan: A5
# REPLACES Yam_RDA.R and panel A of Combined_RDA_GenomicOffset_Figure.R.
# Fixes relative to Yam_RDA.R:
#   - env rows taken from the matched master object (original used unmatched env_data)
#   - VIF pruning instead of |r| > 0.9 filter computed on a different column set
#   - environment conditioned on ancestry (PCs) and geography; varpart with adj. R2
#   - 999 permutations by model, axis and term (original: 100)
# Outputs: Table 2, Fig. 3A (full + partial RDA), within-species pRDA.
# v3: RDA fitted to the full set of genetic PC scores instead of the SNP matrix
#     (mathematically identical results; avoids 'vector memory limit reached').
################################################################################
source("00_config_and_helpers.R"); need("vegan"); library(vegan); M <- load_master(); meta <- M$meta

run_rda <- function(ids, label, ancestry = TRUE) {
  m <- meta[match(ids, meta$SampleID), ]
  Y <- maf_filter(M$G[ids, M$thin, drop = FALSE])       # thinned SNPs (see CFG$thin_bp)
  if (CFG$rda_max_snps > 0 && ncol(Y) > CFG$rda_max_snps) {
    set.seed(CFG$seed); Y <- Y[, sort(sample(ncol(Y), CFG$rda_max_snps)), drop = FALSE] }
  storage.mode(Y) <- "double"; Y[is.na(Y)] <- 0
  n_snps <- ncol(Y)
  # Memory/speed: RDA depends on Y only through its Euclidean geometry, so replacing the
  # (n x SNPs) matrix by its full set of principal-component scores (n x <= n-1) gives
  # IDENTICAL eigenvalues, R2, adjusted R2, variance partitioning, F statistics,
  # permutation P-values and site scores, using a few MB instead of many GB.
  pcY <- prcomp(Y, center = TRUE, scale. = FALSE); keep <- pcY$sdev > 1e-8 * max(pcY$sdev)
  Y <- pcY$x[, keep, drop = FALSE]; rm(pcY); invisible(gc())
  E <- as.data.frame(scale(M$E[ids, M$vif_keep])); vp <- vif_prune(E); E <- E[, vp$keep, drop = FALSE]
  S <- data.frame(lat = scale(m$Latitude)[, 1], lon = scale(m$Longitude)[, 1])
  if (label == "praehensilis") A <- data.frame(lineage = as.numeric(m$Lineage %in% "Cameroonian"))
  else { A <- as.data.frame(Y[, 1:min(CFG$n_pc, ncol(Y) - 1), drop = FALSE]) }   # Y now holds genetic PC scores
  A <- A[, vapply(A, function(v) sd(v) > 0, logical(1)), drop = FALSE]; if (!ncol(A)) A <- NULL
  if (!ancestry) A <- NULL
  D <- if (is.null(A)) cbind(E, S) else cbind(E, S, A)
  cond <- paste(c(names(S), if (!is.null(A)) names(A)), collapse = " + ")
  full <- rda(as.formula(paste("Y ~", paste(names(D), collapse = " + "))), data = D)
  penv <- rda(as.formula(paste("Y ~", paste(names(E), collapse = " + "), "+ Condition(", cond, ")")), data = D)  # env | geo + ancestry
  vpart <- if (!is.null(A)) varpart(Y, E, A, S) else varpart(Y, E, S)
  msg(label, ": permutation tests (", CFG$perms, ")...")
  a_mod  <- anova.cca(penv, permutations = CFG$perms, parallel = CFG$rda_cores)
  a_axis <- anova.cca(penv, by = "axis", permutations = CFG$perms, parallel = CFG$rda_cores)
  a_term <- anova.cca(penv, by = "margin", permutations = CFG$perms, parallel = CFG$rda_cores)
  fr <- vpart$part$indfract
  list(full = full, penv = penv, vpart = vpart, E = E,
       summary = data.frame(set = label, n = length(ids), n_SNPs = n_snps, n_env = ncol(E),
         full_R2 = RsquareAdj(full)$r.squared, full_adjR2 = RsquareAdj(full)$adj.r.squared,
         env_given_geo_anc_adjR2 = RsquareAdj(penv)$adj.r.squared, env_given_geo_anc_P = a_mod$`Pr(>F)`[1],
         pure_env = fr$Adj.R.square[1], pure_anc = if (!is.null(A)) fr$Adj.R.square[2] else NA,
         pure_geo = fr$Adj.R.square[if (!is.null(A)) 3 else 2],
         n_sig_axes = sum(a_axis$`Pr(>F)` < 0.05, na.rm = TRUE),
         sig_terms = paste(rownames(a_term)[which(a_term$`Pr(>F)` < 0.05)], collapse = "; ")),
       axis = data.frame(set = label, a_axis), term = data.frame(set = label, a_term),
       fractions = data.frame(set = label, fraction = rownames(fr), fr))
}
S <- subsets(M)
out <- lapply(names(S), function(s) run_rda(S[[s]], s)); names(out) <- names(S)
write.csv(bind_rows(lapply(out, `[[`, "summary")), tab("Table2_pRDA_summary.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(out, `[[`, "fractions")), tab("Table2_varpart_fractions.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(out, `[[`, "axis")), tab("Table2_tests_by_axis.csv"), row.names = TRUE)
write.csv(bind_rows(lapply(out, `[[`, "term")), tab("Table2_tests_by_term.csv"), row.names = TRUE)
saveRDS(out, res("rda_models.rds")); print(bind_rows(lapply(out, `[[`, "summary")))

# ---- Fig. 3A: full RDA vs partial RDA (environment | ancestry + geography) ----
plot_rda <- function(mod, title) {
  ve <- round(100 * mod$CCA$eig[1:2] / mod$tot.chi, 1)
  si <- as.data.frame(scores(mod, 1:2, "sites")); si$group <- meta$Population.assignment[match(rownames(si), meta$SampleID)]
  bp <- as.data.frame(scores(mod, 1:2, "bp")) * 3; bp$var <- rownames(bp)
  ggplot(si, aes(RDA1, RDA2, colour = group)) + geom_point(alpha = .8) +
    geom_segment(data = bp, aes(0, 0, xend = RDA1, yend = RDA2), inherit.aes = FALSE, colour = "red3",
                 arrow = arrow(length = unit(.15, "cm"))) +
    geom_text(data = bp, aes(RDA1, RDA2, label = var), inherit.aes = FALSE, size = 2.3) +
    labs(title = title, x = paste0("RDA1 (", ve[1], "% of total)"), y = paste0("RDA2 (", ve[2], "% of total)"), colour = "Group") +
    theme_bw()
}
p_full <- plot_rda(out$full$full, "A1  Full RDA (environment + geography + ancestry)")
p_part <- plot_rda(out$full$penv, "A2  Partial RDA: environment | ancestry + geography")
ggsave(fig("Fig3A_full_RDA.png"), p_full, width = 8, height = 6, dpi = 300)
ggsave(fig("Fig3A_partial_RDA.png"), p_part, width = 8, height = 6, dpi = 300)
# NOTE: report axis % consistently (here % of TOTAL variance; the original caption/text mixed
# "% of constrained variance"). State which one in the caption.
