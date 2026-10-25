################################################################################
# 07_EGWAS_revised.R                                       Analysis Plan: A6a, A6b
# REPLACES Yam_egwas_rrblup_V2.R, Yam_GWAS_rrBLUP.R and Yam_EGWAS.R (GAPIT), and
# the separate GEMMA run on ICs.
#   - one mixed-model scan (EMMAX/P3D, as rrBLUP::GWAS with P3D = TRUE) returning
#     beta, SE, allele frequency and PVE for every SNP
#   - Q + K model (kinship + CFG$n_pc PCs) for the full panel
#   - within-species scans for D. rotundata and D. praehensilis (lineage covariate)
#   - latitude/longitude are NOT environmental variables (excluded in 01)
#   - FDR = CFG$fdr (Benjamini-Hochberg); genomic inflation (lambda) per trait and set
#   - ICs (icaimax, 3 components on the 19 bioclim variables) computed and saved here
# v2: covariate matrix made full rank in emmax_scan (constant lineage covariate no
#     longer stops the run); sets smaller than CFG$min_n_within are skipped with a
#     report; NA p-values (monomorphic SNPs) removed before FDR; Fig. S6 with
#     chromosome axis, any number of contigs, and QQ plots.
# v3: responses rank-inverse-normal transformed (CFG$int_transform); REML bounded
#     (CFG$min_ve_ratio); PVE = t^2/(t^2+df) (bounded); 'boundary' flag in Table S7.
# Outputs: Table S5 (IC loadings), Table S18 (+ S4b lambda/counts), EGWAS summary,
#          per-SNP p-values for WZA (08), Fig. S6.
################################################################################
RUN_SCAN <- TRUE    # set FALSE to redraw Fig. S6 from saved p-values only
source("00_config_and_helpers.R"); need(c("rrBLUP", "ica")); M <- load_master(); meta <- M$meta

# ---- traits: all env variables + IC1-3 ---------------------------------------
ic_file <- res("ICs.rds")
if (RUN_SCAN || !file.exists(ic_file)) {
  bio <- grep("^bio", colnames(M$E), ignore.case = TRUE, value = TRUE)
  ic <- ica::icaimax(scale(M$E[, bio]), nc = 3, center = TRUE)
  IC <- ic$S; colnames(IC) <- paste0("IC", 1:3); rownames(IC) <- rownames(M$E)
  write.csv(round(cor(M$E[, bio], IC), 3), tab("TableS5_IC_loadings.csv"))
  saveRDS(IC, ic_file)
} else IC <- readRDS(ic_file)
Y <- cbind(M$E, IC)

# ---- sample sets ------------------------------------------------------------------
S <- subsets(M)
msg("Set sizes:", paste(names(S), lengths(S), sep = " = ", collapse = "; "))
small <- names(S)[lengths(S) < CFG$min_n_within]
if (length(small)) {
  msg("WARNING: skipping set(s) with fewer than", CFG$min_n_within, "samples:", paste(small, collapse = ", "))
  print(table(meta$Species_group, meta$NonAdmixed, dnn = c("Species_group", "NonAdmixed"), useNA = "ifany"))
  S <- S[setdiff(names(S), small)]
}

# ---- scans ----------------------------------------------------------------------
if (RUN_SCAN) {
all_hits <- list(); lam <- list()
for (s in names(S)) {
  ids <- intersect(S[[s]], rownames(M$G)); m <- meta[match(ids, meta$SampleID), ]
  G <- maf_filter(M$G[ids, ]); K <- kinship(M$G[ids, ])
  X <- switch(s,
    full         = M$PCs[ids, 1:CFG$n_pc],
    praehensilis = { lin <- as.numeric(m$Lineage %in% "Cameroonian")
                     msg("praehensilis lineages: Cameroonian =", sum(lin), "; other =", sum(1 - lin))
                     matrix(lin, ncol = 1, dimnames = list(NULL, "Cameroonian")) },
    rotundata    = prcomp(G)$x[, 1:2])          # 2 within-species PCs; set NULL for K-only
  X <- full_rank_X(X, length(ids))[, -1, drop = FALSE]; if (!ncol(X)) X <- NULL   # reduce once per set
  msg("Set", s, ":", length(ids), "samples,", ncol(G), "SNPs,", if (is.null(X)) 0 else ncol(X), "covariates")
  P <- list()
  for (tr in colnames(Y)) {
    y <- if (isTRUE(CFG$int_transform)) rint(Y[ids, tr]) else as.numeric(scale(Y[ids, tr]))
    if (!is.finite(sd(Y[ids, tr])) || sd(Y[ids, tr]) == 0) { msg("  skip", tr, "(no variation in", s, ")"); next }
    r <- emmax_scan(y, G, K, X)
    r$q <- p.adjust(r$p, "BH"); r$trait <- tr; r$set <- s          # p.adjust ignores NA p-values
    lam[[paste(s, tr)]] <- data.frame(set = s, trait = tr, lambda = lambda_gc(r$p), h2_null = r$h2_null[1],
                                      n_tested = sum(is.finite(r$p)), n_sig = sum(r$q < CFG$fdr, na.rm = TRUE),
                                      boundary = r$h2_null[1] > 0.94)
    P[[tr]] <- r[, c("SNP", "beta", "se", "p", "af", "PVE")]
    all_hits[[paste(s, tr)]] <- r %>% filter(!is.na(q), q < CFG$fdr)
  }
  saveRDS(P, res(paste0("egwas_pvalues_", s, ".rds")))
}
hits <- bind_rows(all_hits)
if (nrow(hits)) {
  hits <- hits %>% left_join(M$snp[, c("SNP", "CHROM", "POS")], by = "SNP")
  hits$beta_original_units <- if (isTRUE(CFG$int_transform)) NA_real_ else hits$beta * apply(Y, 2, sd)[hits$trait]   # beta is in INT/SD units
} else hits <- data.frame(SNP = character(), beta = numeric(), se = numeric(), p = numeric(), af = numeric(),
                          PVE = numeric(), q = numeric(), trait = character(), set = character(),
                          CHROM = character(), POS = numeric(), beta_original_units = numeric())
lam <- bind_rows(lam)
write.csv(hits, tab("TableS18_EGWAS_significant_with_effects.csv"), row.names = FALSE)
write.csv(lam, tab("TableS7_lambda_and_counts.csv"), row.names = FALSE)

# ---- summaries used in the Results ------------------------------------------
summ <- hits %>% group_by(set) %>% summarise(associations = n(), unique_markers = n_distinct(SNP),
          median_PVE_pct = 100 * median(PVE), min_PVE_pct = 100 * min(PVE), max_PVE_pct = 100 * max(PVE), .groups = "drop")
full_m <- unique(hits$SNP[hits$set == "full"])
summ$recovered_from_full <- sapply(summ$set, function(s) if (s == "full") NA else length(intersect(full_m, hits$SNP[hits$set == s])))
summ$lambda_range <- sapply(summ$set, function(s) paste(round(range(lam$lambda[lam$set == s], na.rm = TRUE), 2), collapse = "-"))
write.csv(summ, tab("EGWAS_summary_by_set.csv"), row.names = FALSE); print(summ)
if (length(full_m)) {
  pct_multi <- hits %>% filter(set == "full") %>% count(SNP) %>% summarise(p = 100 * mean(n > 1)) %>% pull(p)
  msg("Full panel: % markers associated with >1 variable =", round(pct_multi, 1))
}
}

# ---- Fig. S6: Manhattan + QQ for IC1-3 (full panel) --------------------------
pf <- res("egwas_pvalues_full.rds")
if (file.exists(pf)) {
  P <- readRDS(pf)
  for (tr in paste0("IC", 1:3)) {
    if (is.null(P[[tr]])) { msg("No results for", tr); next }
    d <- P[[tr]] %>% filter(is.finite(p), p > 0) %>%
      left_join(M$snp[, c("SNP", "CHROM", "POS")], by = "SNP") %>%
      arrange(CHROM, POS) %>% mutate(idx = row_number(), q = p.adjust(p, "BH"))
    sig <- d$q < CFG$fdr; thr <- if (any(sig)) max(d$p[sig]) else NA; lam_tr <- lambda_gc(d$p)
    nchr <- n_distinct(d$CHROM)
    ax <- d %>% group_by(CHROM) %>% summarise(mid = (min(idx) + max(idx)) / 2, .groups = "drop")
    pm <- ggplot(d, aes(idx, -log10(p), colour = CHROM)) + geom_point(size = .4, show.legend = FALSE) +
      scale_colour_manual(values = rep(c("black", "grey60"), length.out = nchr)) +
      scale_x_continuous(breaks = ax$mid, labels = seq_len(nrow(ax))) +
      { if (!is.na(thr)) geom_hline(yintercept = -log10(thr), colour = "red", linetype = 2) } +
      labs(title = NULL,
           x = "Chromosome", y = expression(-log[10](p))) + theme_bw() + theme(axis.text.x = element_text(size = 6))
    qq <- data.frame(exp = -log10(ppoints(nrow(d))), obs = -log10(sort(d$p)))
    pq <- ggplot(qq, aes(exp, obs)) + geom_point(size = .5) + geom_abline(slope = 1, intercept = 0, colour = "red") +
      labs(title = NULL, x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p))) + theme_bw()
    save_panel(pm, paste0("FigS6_", tr, "_manhattan")); save_panel(pq, paste0("FigS6_", tr, "_qq"))
    if (requireNamespace("patchwork", quietly = TRUE)) {
      ggsave(fig(paste0("FigS6_", tr, ".png")), patchwork::wrap_plots(pm, pq, widths = c(3, 1)) + patchwork::plot_annotation(tag_levels = "A"), width = 14, height = 3.8, dpi = 300)
    } else {
      ggsave(fig(paste0("FigS6_manhattan_", tr, ".png")), pm, width = 11, height = 3.5, dpi = 300)
      ggsave(fig(paste0("FigS6_qq_", tr, ".png")), pq, width = 4, height = 4, dpi = 300)
    }
  }
  write.csv(data.frame(chromosome_label = seq_len(nrow(ax)), contig = ax$CHROM), tab("FigS6_chromosome_labels.csv"), row.names = FALSE)
} else msg("No full-panel p-values found - run with RUN_SCAN <- TRUE")
