################################################################################
# 09_LFMM.R                                                Analysis Plan: A6e
# LFMM2 (LEA) association results - described in the original Methods but never
# reported. Full panel K = 7; within species K from CFG$K_lfmm_within (check with
# the snmf cross-entropy block below). Also saves the LFMM input used by 11.
# Output: Table S8.  Adapted from the lfmm2 block of
# General_genetic_offset_script_with_regional_plots.R.
################################################################################
source("00_config_and_helpers.R"); need("LEA"); M <- load_master()
S <- subsets(M); hits <- read.csv(tab("TableS18_EGWAS_significant_with_effects.csv")); out <- list()
for (s in names(S)) {
  ids <- S[[s]]; G <- maf_filter(M$G[ids, ]); E <- scale(M$E[ids, M$vif_keep])
  f <- res(paste0("geno_", s, ".lfmm")); LEA::write.lfmm(G, f)
  K <- if (s == "full") CFG$K_lfmm_full else CFG$K_lfmm_within[[s]]
  # --- optional check of K (uncomment):
  # sn <- LEA::snmf(LEA::lfmm2geno(f), K = 1:8, entropy = TRUE, repetitions = 3, project = "new"); plot(sn)
  mod <- LEA::lfmm2(input = f, env = E, K = K, effect.sizes = TRUE)
  pv <- LEA::lfmm2.test(mod, input = f, env = E, full = FALSE, genomic.control = TRUE)
  pm <- pv$pvalues; if (is.null(dim(pm))) pm <- matrix(pm, ncol = 1)
  if (nrow(pm) != ncol(G)) pm <- base::t(pm)                       # orientation differs among LEA versions
  colnames(pm) <- colnames(E)
  q <- apply(pm, 2, p.adjust, method = "BH")
  sig <- which(q < CFG$fdr, arr.ind = TRUE)
  tabs <- data.frame(set = s, K = K, SNP = colnames(G)[sig[, 1]], trait = colnames(E)[sig[, 2]], p = pm[sig], q = q[sig])
  tabs$shared_with_EGWAS <- paste(tabs$SNP, tabs$trait) %in% paste(hits$SNP[hits$set == s], hits$trait[hits$set == s])
  out[[s]] <- tabs; saveRDS(list(mod = mod, snps = colnames(G), env = colnames(E)), res(paste0("lfmm_", s, ".rds")))
}
s10 <- bind_rows(out); write.csv(s10, tab("TableS8_LFMM_significant.csv"), row.names = FALSE)
print(s10 %>% group_by(set) %>% summarise(n_SNP_trait = n(), n_SNPs = n_distinct(SNP), shared = sum(shared_with_EGWAS)))
