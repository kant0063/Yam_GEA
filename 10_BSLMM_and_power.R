################################################################################
# 10_BSLMM_and_power.R                                     Analysis Plan: A7
# (a) BSLMM in GEMMA (Zhou et al. 2013): PVE, PGE, n_gamma for IC1-3 and key
#     variables -> Table S11. Needs the GEMMA binary (CFG$gemma_bin).
# (b) Power simulation on the real genotypes/kinship -> Fig. S18.
# New script. BSLMM: 3 chains x 10 M steps per trait; run on a server.
#
# v4: power simulation robust to REML failures (see section b).
# v3: GEMMA path checks; robust .hyp.txt reader (trailing tabs); summary no longer
# overwrites n_gamma before its quantiles; per-chain table; switches below so the
# summary can be rebuilt without rerunning GEMMA or the (slow) power simulation.
################################################################################
RUN_BSLMM <- TRUE    # GEMMA runs are skipped automatically when their .hyp.txt exists
RUN_POWER <- TRUE    # set FALSE to rebuild Table S11 only
source("00_config_and_helpers.R"); M <- load_master(); meta <- M$meta; S <- subsets(M)
msg("Set sizes:", paste(names(S), lengths(S), sep = " = ", collapse = "; "))
small <- names(S)[lengths(S) < CFG$min_n_within]
if (length(small)) {
  msg("WARNING: skipping set(s) with fewer than", CFG$min_n_within, "samples:", paste(small, collapse = ", "),
      "- check Population.assignment / Species_group / NonAdmixed in the master object:")
  print(table(meta$Species_group, meta$NonAdmixed, dnn = c("Species_group", "NonAdmixed"), useNA = "ifany"))
  S <- S[setdiff(names(S), small)]
}
IC <- readRDS(res("ICs.rds")); Y <- cbind(M$E, IC)
traits <- na.omit(sapply(CFG$heavy_traits, function(t) find_var(colnames(Y), t)))
dir.create(res("gemma"), showWarnings = FALSE)

# ---- GEMMA executable: absolute path, executable, and actually runs -------------
# (system() runs from a shell that does not see your PATH, and the command below
#  does "cd results/<set>/gemma", so a bare "gemma" or a relative path fails.)
if (RUN_BSLMM) {
gemma <- normalizePath(path.expand(CFG$gemma_bin), mustWork = FALSE)
if (!file.exists(gemma)) stop("GEMMA not found at ", gemma, " - set CFG$gemma_bin in 00 or Sys.setenv(YAM_GEMMA = ...)")
if (file.access(gemma, 1) != 0) { msg("Making GEMMA executable"); Sys.chmod(gemma, "755") }
test <- suppressWarnings(system2(gemma, "-license", stdout = TRUE, stderr = TRUE))
if (!is.null(attr(test, "status")) && attr(test, "status") != 0 && !any(grepl("GEMMA", test)))
  stop("GEMMA exists but will not run:\n", paste(test, collapse = "\n"),
       "\nOn macOS try in Terminal:  xattr -d com.apple.quarantine ", gemma,
       "\nIf it is a Linux binary, install a macOS build instead (e.g. conda install -c bioconda gemma).")
msg("Using GEMMA:", gemma)
}

# ---- (a) BSLMM -------------------------------------------------------------
if (RUN_BSLMM) {
write_bimbam <- function(G, snp, file) {     # mean-genotype format: id, alt, ref, dosages
  s <- snp[match(colnames(G), snp$SNP), ]
  fwrite(data.table(s$SNP, s$ALT, s$REF, base::t(G)), file, col.names = FALSE, sep = ",")
}
# GEMMA .hyp.txt: whitespace-separated, often with a trailing tab on every line.
# Read it by splitting on whitespace so columns cannot shift, and check the header.
read_hyp <- function(f) {
  L <- readLines(f); L <- L[nzchar(trimws(L))]
  hdr <- strsplit(trimws(L[1]), "[[:space:]]+")[[1]]
  vals <- lapply(strsplit(trimws(L[-1]), "[[:space:]]+"), as.numeric)
  ok <- lengths(vals) == length(hdr)
  if (!all(ok)) warning(f, ": ", sum(!ok), " malformed rows dropped")
  d <- as.data.frame(do.call(rbind, vals[ok])); names(d) <- hdr
  need_cols <- c("pve", "pge", "n_gamma")
  if (!all(need_cols %in% names(d))) stop(f, ": expected columns ", paste(need_cols, collapse = ", "),
                                           "; found ", paste(hdr, collapse = ", "))
  d
}
hyp <- list()
for (s in intersect(c("full", "rotundata"), names(S))) {
  ids <- S[[s]]; gf <- res(paste0("gemma/geno_", s, ".txt"))
  if (!file.exists(gf)) { G <- maf_filter(M$G[ids, ]); write_bimbam(G, M$snp, gf) }
  for (tr in traits) {
    pf <- res(paste0("gemma/pheno_", s, "_", tr, ".txt")); fwrite(data.table(scale(Y[ids, tr])), pf, col.names = FALSE)
    for (ch in 1:3) {
      o <- paste0("bslmm_", s, "_", tr, "_c", ch)
      hyp_file <- res(paste0("gemma/output/", o, ".hyp.txt"))
      cmd <- sprintf("cd %s && %s -g %s -p %s -bslmm 1 -w 1000000 -s 10000000 -seed %d -o %s",
                     shQuote(res("gemma")), shQuote(gemma), shQuote(basename(gf)), shQuote(basename(pf)), CFG$seed + ch, o)
      if (!file.exists(hyp_file)) {
        msg(cmd); rc <- system(cmd)
        if (rc != 0 || !file.exists(hyp_file))
          stop("GEMMA failed for ", o, " (exit ", rc, "). See ", res(paste0("gemma/output/", o, ".log.txt")))
      }
      h <- read_hyp(hyp_file)
      hyp[[o]] <- data.frame(set = s, trait = tr, chain = ch, pve = h$pve, pge = h$pge, n_gamma = h$n_gamma)
    }
  }
}
H <- bind_rows(hyp)
msg("Posterior samples read:", nrow(H), "; NA values: pve", sum(is.na(H$pve)), "pge", sum(is.na(H$pge)), "n_gamma", sum(is.na(H$n_gamma)))
# quantiles are computed BEFORE the medians so no summary overwrites a column it still needs
s11 <- H %>% group_by(set, trait) %>%
  summarise(PVE_lo = quantile(pve, .025, na.rm = TRUE), PVE_hi = quantile(pve, .975, na.rm = TRUE),
            PGE_lo = quantile(pge, .025, na.rm = TRUE), PGE_hi = quantile(pge, .975, na.rm = TRUE),
            n_gamma_lo = quantile(n_gamma, .025, na.rm = TRUE), n_gamma_hi = quantile(n_gamma, .975, na.rm = TRUE),
            PVE = median(pve, na.rm = TRUE), PGE = median(pge, na.rm = TRUE), n_gamma_median = median(n_gamma, na.rm = TRUE),
            n_samples = n(), .groups = "drop") %>%
  dplyr::select(set, trait, PVE, PVE_lo, PVE_hi, PGE, PGE_lo, PGE_hi, n_gamma_median, n_gamma_lo, n_gamma_hi, n_samples)
# chain agreement: medians per chain should be similar
chains <- H %>% group_by(set, trait, chain) %>% summarise(pve = median(pve, na.rm = TRUE), pge = median(pge, na.rm = TRUE),
                                                         n_gamma = median(n_gamma, na.rm = TRUE), .groups = "drop")
write.csv(chains, tab("TableS11b_BSLMM_by_chain.csv"), row.names = FALSE); print(chains)
write.csv(s11, tab("TableS11_BSLMM.csv"), row.names = FALSE); print(s11)
# Check chain agreement: medians of pve/pge across the three chains should be close for each trait.
}

# ---- (b) Power simulation ---------------------------------------------------
# v4: robust to occasional REML failures ("X not full rank" from rrBLUP when a
# simulated response is not finite or the kinship is not positive definite):
#   - kinship made positive definite (eigen-decomposition, small eigenvalues floored)
#   - causal SNP drawn only from SNPs that vary within the set
#   - each replicate is checked and wrapped in tryCatch; if REML fails, the scan uses
#     the true simulated variance components (P3D) instead, and failures are counted
if (RUN_POWER) {
pd_root <- function(K) {                        # returns L with L %*% base::t(L) ~ K (PSD fix)
  e <- eigen((K + base::t(K)) / 2, symmetric = TRUE); ev <- pmax(e$values, 1e-6 * max(e$values))
  list(L = e$vectors %*% diag(sqrt(ev)), K = e$vectors %*% diag(ev) %*% base::t(e$vectors))
}
h2_bg <- 0.3
pw <- list(); fails <- list()
for (s in names(S)) {
  ids <- intersect(S[[s]], rownames(M$G)); n <- length(ids)
  G <- maf_filter(M$G[ids, ]); storage.mode(G) <- "double"; G[is.na(G)] <- 0
  sdj <- apply(G, 2, sd); cand <- which(is.finite(sdj) & sdj > 0)
  pk <- pd_root(kinship(M$G[ids, ])); K <- pk$K; L <- pk$L; dimnames(K) <- list(ids, ids)
  msg("Power:", s, "| n =", n, "| SNPs =", ncol(G), "| candidate causal SNPs =", length(cand))
  for (pve in CFG$power_pve) {
    n_fail <- 0
    det <- vapply(seq_len(CFG$power_reps), function(i) {
      j <- cand[sample.int(length(cand), 1)]
      x <- G[, j] - mean(G[, j])
      g <- as.numeric(L %*% rnorm(n)); g <- g / sd(g) * sqrt(h2_bg)
      e <- rnorm(n); e <- e / sd(e) * sqrt(1 - h2_bg - pve)
      y <- x / sd(x) * sqrt(pve) + g + e
      if (!all(is.finite(y))) { n_fail <<- n_fail + 1; return(NA) }
      r <- tryCatch(emmax_scan(y, G, K),
                    error = function(err) { n_fail <<- n_fail + 1
                      emmax_scan(y, G, K, fit = list(Vu = h2_bg / mean(diag(K)), Ve = 1 - h2_bg - pve)) })
      isTRUE(p.adjust(r$p, "BH")[j] < CFG$fdr)
    }, logical(1))
    pw[[paste(s, pve)]] <- data.frame(set = s, n = n, PVE = pve, power = mean(det, na.rm = TRUE),
                                      reps = sum(!is.na(det)), REML_fallbacks = n_fail)
    msg(s, "PVE", pve, "power", round(mean(det, na.rm = TRUE), 3), if (n_fail) paste("(REML fallbacks:", n_fail, ")") else "")
  }
}
pw <- bind_rows(pw); write.csv(pw, tab("FigS18_power.csv"), row.names = FALSE); print(pw)
ggsave(fig("FigS18_power.png"), ggplot(pw, aes(100 * PVE, power, colour = set)) + geom_line() + geom_point() +
  geom_hline(yintercept = .8, linetype = 2) + labs(x = "Variance explained by the locus (%)",
  y = paste0("Detection rate (FDR ", CFG$fdr, ")"), colour = NULL) + theme_bw(), width = 6, height = 4, dpi = 300)
}
