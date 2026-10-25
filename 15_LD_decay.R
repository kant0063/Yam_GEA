################################################################################
# 15_LD_decay.R                                            Fig. S12; window sizes
# No original script was uploaded for Fig. S12. Recomputes LD decay (r^2 between SNP
# pairs <= CFG$ld_max_dist apart) for the full panel and within D. rotundata, fits
# r2 = c + a*exp(-b*d), and reports the half-decay distance (d where r2 falls halfway
# from its initial value to background) and the distance to r2 = 0.2 (NA if r2 never
# reaches 0.2 - this resolves the Methods/caption inconsistency).
# Outputs: FigS12_LD_decay.png, LD_decay_summary.csv (used for gene / WZA windows).
################################################################################
source("00_config_and_helpers.R")
# LD needs dense markers: use the relaxed (imputed, 1.57 M SNP) master if it exists, else the strict one
rel <- file.path("results", "relaxed", "master.rds")
M <- if (file.exists(rel)) { msg("Using relaxed SNP set for LD decay"); readRDS(rel) } else load_master()
S <- subsets(M)
ld_pairs <- function(G, snp) {
  snp <- snp[match(colnames(G), snp$SNP), ]; out <- list()
  anchors <- sample(ncol(G), min(CFG$ld_anchors, ncol(G)))
  for (a in anchors) {
    j <- which(snp$CHROM == snp$CHROM[a] & snp$POS > snp$POS[a] & snp$POS - snp$POS[a] <= CFG$ld_max_dist)
    if (length(j)) out[[length(out) + 1]] <- data.frame(d = snp$POS[j] - snp$POS[a],
                                                       r2 = as.numeric(suppressWarnings(cor(G[, a], G[, j, drop = FALSE])))^2)
  }
  bind_rows(out) %>% filter(is.finite(r2))
}
summ <- list(); plots <- list()
for (s in c("full", "rotundata")) {
  Gs <- maf_filter(M$G[intersect(S[[s]], rownames(M$G)), ]); storage.mode(Gs) <- "double"; Gs[is.na(Gs)] <- 0
  P <- ld_pairs(Gs, M$snp)
  msg(s, ":", ncol(Gs), "SNPs;", nrow(P), "SNP pairs within", CFG$ld_max_dist / 1000, "kb from", CFG$ld_anchors, "anchors")
  if (nrow(P) < 1000) warning("Few SNP pairs (", nrow(P), ") - LD decay estimate unreliable; check SNP positions")
  P$bin <- cut(P$d, unique(c(0, exp(seq(log(100), log(CFG$ld_max_dist), length = 60)))))
  B <- P %>% group_by(bin) %>% summarise(d = median(d), r2 = mean(r2), .groups = "drop")
  fit <- tryCatch(nls(r2 ~ c0 + a * exp(-b * d), data = B, start = list(c0 = min(B$r2), a = max(B$r2) - min(B$r2), b = 1e-4)),
                  error = function(e) NULL)
  if (!is.null(fit)) { cf <- coef(fit); half <- log(2) / cf["b"]
    d02 <- if (cf["c0"] + cf["a"] > 0.2 && cf["c0"] < 0.2) -log((0.2 - cf["c0"]) / cf["a"]) / cf["b"] else NA
  } else { cf <- c(c0 = NA, a = NA, b = NA); half <- NA; d02 <- NA }
  summ[[s]] <- data.frame(set = s, n = length(S[[s]]), pairs = nrow(P), r2_initial = cf["c0"] + cf["a"], r2_background = cf["c0"],
                          half_decay_bp = half, dist_to_r2_0.2_bp = d02)
  B$set <- s; B$fit <- if (!is.null(fit)) predict(fit, B) else NA_real_; plots[[s]] <- B
}
summ <- bind_rows(summ); write.csv(summ, tab("LD_decay_summary.csv"), row.names = FALSE); print(summ)
pp <- ggplot(bind_rows(plots), aes(d / 1000, r2, colour = set)) + geom_point(size = 1) + geom_line(aes(y = fit), na.rm = TRUE) +
  scale_x_log10() + labs(x = "Distance (kb, log scale)", y = expression(Mean~r^2)) + theme_bw()
ggsave(fig("FigS12_LD_decay.png"), pp, width = 6, height = 4, dpi = 300)
