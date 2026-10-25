################################################################################
# 13_genomic_prediction_revised.R                          Analysis Plan: A9
# REPLACES xval_kfold_functions.R (and the GEAV prediction step that used the
# 50-accession core as training set).
# Problems in the original functions:
#   - CV folds were drawn only within the 50-line core (y.trainset); accuracy was
#     computed on the held-out core lines PLUS all 117 non-core lines pooled.
#   - k.xval.EXP called kin.blup(..., GAUSS = TRUE), i.e. the "exponential" model was
#     fitted as a Gaussian kernel, so two of the three models may be the same model.
#   - no baseline: environment of origin is largely predictable from ancestry, so
#     random-CV accuracy and rankings mostly reproduce species/region (circularity).
# v3: the core is used again as a training set (core -> rest, as in the EGS framework),
#     benchmarked against random and group-stratified training sets of the same size,
#     with the structure-only baseline; Table S14 built from genotypes only if no core file.
# v2: mixed.solve beta is a 1x1 matrix -> coerced (fixes 'non-conformable arrays');
#     drop = FALSE for single-accession folds; prediction length checked.
# This script: all 167 accessions; random k-fold, leave-one-group-out and spatially
# blocked CV; structure-only baseline; rrBLUP, GBLUP-Gaussian, GBLUP-exponential,
# BayesB, RKHS, random forest (+ optional GenoBridge hook); within-group GEAVs.
# Outputs: Fig. S9, Fig. S10, Table S16 (GEAV columns), candidate parents.
################################################################################
source("00_config_and_helpers.R"); need(c("rrBLUP", "BGLR", "ranger")); M <- load_master(); meta <- M$meta
IC <- readRDS(res("ICs.rds")); Y <- scale(cbind(M$E, IC)); ids <- rownames(Y)
G <- maf_filter(M$G[ids, ]); Gc <- scale(G, scale = FALSE); PC <- M$PCs[ids, 1:CFG$n_pc]
D <- as.matrix(dist(Gc)); D <- D / max(D); theta <- median(D[upper.tri(D)])
K_gau <- exp(-(D / theta)^2); K_exp <- exp(-D / theta)
heavy <- na.omit(sapply(CFG$heavy_traits, function(t) find_var(colnames(Y), t)))

# ---- fold definitions ------------------------------------------------------
folds_random <- function() { f <- sample(rep(1:CFG$cv_k, length.out = length(ids))); split(seq_along(ids), f) }
folds_logo   <- split(seq_along(ids), meta$Population.assignment[match(ids, meta$SampleID)])
km <- kmeans(meta[match(ids, meta$SampleID), c("Longitude", "Latitude")], CFG$n_spatial_blocks, nstart = 50)
folds_space  <- split(seq_along(ids), km$cluster)
# (blockCV::cv_spatial() is an alternative for spatial blocks; k-means on coordinates keeps it dependency-free)

# ---- models: each returns predictions for test indices -------------------------
fit_pred <- function(model, y, tr, te) {
  # rrBLUP::mixed.solve returns beta as a 1 x 1 matrix: always use b0 <- as.numeric(beta)[1].
  # drop = FALSE everywhere so a fold with a single test accession still works.
  yNA <- y; yNA[te] <- NA
  pred <- switch(model,
    baseline = { d_tr <- data.frame(y = y[tr], PC[tr, , drop = FALSE]); d_te <- data.frame(PC[te, , drop = FALSE])
                 predict(lm(y ~ ., data = d_tr), newdata = d_te) },
    rrBLUP   = { f <- rrBLUP::mixed.solve(y[tr], Z = Gc[tr, , drop = FALSE])
                 as.numeric(f$beta)[1] + as.numeric(Gc[te, , drop = FALSE] %*% f$u) },
    GBLUP_gauss = { f <- rrBLUP::mixed.solve(yNA, K = K_gau); as.numeric(f$beta)[1] + as.numeric(f$u)[te] },
    GBLUP_exp   = { f <- rrBLUP::mixed.solve(yNA, K = K_exp); as.numeric(f$beta)[1] + as.numeric(f$u)[te] },
    BayesB = { f <- BGLR::BGLR(yNA, ETA = list(list(X = Gc, model = "BayesB")), nIter = CFG$bglr_iter,
                               burnIn = CFG$bglr_burn, verbose = FALSE, saveAt = tempfile()); f$yHat[te] },
    RKHS   = { f <- BGLR::BGLR(yNA, ETA = list(list(K = K_gau, model = "RKHS")), nIter = CFG$bglr_iter,
                               burnIn = CFG$bglr_burn, verbose = FALSE, saveAt = tempfile()); f$yHat[te] },
    RF     = { f <- ranger::ranger(x = G[tr, , drop = FALSE], y = y[tr], num.trees = 500, seed = CFG$seed)
               predict(f, G[te, , drop = FALSE])$predictions },
    # GenoBridge = { ... call GenoBridge here with the same tr/te indices if it is included ... },
    stop("unknown model"))
  pred <- as.numeric(pred)
  if (length(pred) != length(te)) stop(model, ": returned ", length(pred), " predictions for ", length(te), " test accessions")
  pred
}
run_cv <- function(model, trait, folds) {
  y <- Y[, trait]; pred <- rep(NA_real_, length(y))
  for (te in folds) { tr <- setdiff(seq_along(y), te); pred[te] <- fit_pred(model, y, tr, te) }
  pred
}
grp <- meta$Population.assignment[match(ids, meta$SampleID)]
acc_stats <- function(pred, y) c(r_pooled = cor(pred, y),
  r_within_group = mean(sapply(split(seq_along(y), grp), function(i) if (length(i) > 4) cor(pred[i], y[i]) else NA), na.rm = TRUE))

models_light <- c("baseline", "rrBLUP", "GBLUP_gauss", "GBLUP_exp"); models_heavy <- c("BayesB", "RKHS", "RF")
res_acc <- list(); oof <- list()
for (trait in colnames(Y)) {
  mods <- c(models_light, if (trait %in% heavy) models_heavy)
  for (mo in mods) {
    msg(trait, mo)
    reps <- if (mo %in% models_heavy) 2 else CFG$cv_reps
    pr <- sapply(1:reps, function(i) run_cv(mo, trait, folds_random()))
    res_acc[[paste(trait, mo, "random")]] <- data.frame(trait, model = mo, scheme = "random 10-fold",
      base::t(colMeans(base::t(apply(pr, 2, acc_stats, y = Y[, trait])))))
    oof[[paste(trait, mo)]] <- rowMeans(pr)
    for (sch in c("leave-one-group-out", "spatial blocks")) {
      p <- run_cv(mo, trait, if (sch == "spatial blocks") folds_space else folds_logo)
      res_acc[[paste(trait, mo, sch)]] <- data.frame(trait, model = mo, scheme = sch, base::t(acc_stats(p, Y[, trait])))
    }
  }
}
acc <- bind_rows(res_acc)
base <- acc %>% filter(model == "baseline") %>% dplyr::select(trait, scheme, r_base = r_pooled)
acc <- acc %>% left_join(base, by = c("trait", "scheme")) %>% mutate(gain_over_baseline = r_pooled - r_base)
write.csv(acc, tab("FigS9_prediction_accuracy.csv"), row.names = FALSE)
print(acc %>% group_by(scheme, model) %>% summarise(r = mean(r_pooled), gain = mean(gain_over_baseline), .groups = "drop"))

acc$scheme <- factor(acc$scheme, levels = c("random 10-fold", "leave-one-group-out", "spatial blocks"))
pS14 <- ggplot(acc, aes(trait, r_pooled, colour = model)) + geom_point(position = position_dodge(.6), size = 1.3) +
  facet_wrap(~scheme, ncol = 1) + geom_hline(yintercept = 0, linetype = 2) +
  labs(x = NULL, y = "Prediction accuracy (r)") + theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6))
ggsave(fig("FigS9_prediction_accuracy.png"), pS14, width = 12, height = 9, dpi = 300)

# ---- Core collection: see 13b_core_comparison.R -------------------------------------
# The core-training analysis (Table S14, Table S15, Fig. S9D) is run ONLY by script 13b, which
# compares the original and genotype-based cores and keeps the more accurate one. It is not
# repeated here, so rerunning 13 can no longer overwrite the 13b outputs.

# ---- GEAVs: out-of-fold (random CV), best non-baseline light model per trait,
#      centred on genetic-group means (within-group GEAV) -----------------------
best <- acc %>% filter(scheme == "spatial blocks", model %in% c("rrBLUP", "GBLUP_gauss", "GBLUP_exp")) %>%
  group_by(trait) %>% slice_max(r_pooled, n = 1, with_ties = FALSE)
geav <- sapply(colnames(Y), function(t) oof[[paste(t, best$model[best$trait == t])]])
rownames(geav) <- ids
geav_within <- apply(geav, 2, function(v) v - ave(v, grp))
out <- data.frame(SampleID = ids, group = grp, geav_within, check.names = FALSE)
write.csv(out, tab("TableS16_GEAV_within_group.csv"), row.names = FALSE)
v_ph <- find_var(colnames(Y), VAR_PH); v_pr <- find_var(colnames(Y), VAR_PRECIP)
cand <- list(low_pH = out %>% arrange(.data[[v_ph]]) %>% head(5) %>% dplyr::select(SampleID, group, value = all_of(v_ph)),
             high_BIO12 = out %>% arrange(desc(.data[[v_pr]])) %>% head(5) %>% dplyr::select(SampleID, group, value = all_of(v_pr)))
print(cand); write.csv(bind_rows(cand, .id = "target"), tab("TableS17_candidate_parents_within_group.csv"), row.names = FALSE)
# compare with the original picks: P4937, P4936, P4919, P4921, P4920 (low pH); A3009, A5059, A5061, P599, P424 (BIO12)

pS15 <- plot_geav_heatmap(out, M, cand_ids = unique(unlist(lapply(cand, function(x) x$SampleID))))
write.csv(attr(pS15, "subclusters"), tab("TableS16b_GEAV_subclusters.csv"), row.names = FALSE)
ggsave(fig("FigS10_GEAV_within_group.png"), pS15, width = 10, height = 16, dpi = 300, bg = "white")

