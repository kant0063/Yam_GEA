################################################################################
# 08_WZA.R                                                 Analysis Plan: A6d
# Weighted-Z Analysis (Booker et al. 2024) on the within-species per-SNP EGWAS
# p-values from 07. New script (WZA was introduced but not performed originally).
# Output: Table S9 (top windows per set/trait, overlap with EGWAS hits).
################################################################################
source("00_config_and_helpers.R"); M <- load_master()
hits <- read.csv(tab("TableS4_EGWAS_significant_with_effects.csv"))
out <- list()
for (s in c("rotundata", "praehensilis")) {
  P <- readRDS(res(paste0("egwas_pvalues_", s, ".rds")))
  for (tr in names(P)) {
    d <- P[[tr]] %>% left_join(M$snp[, c("SNP", "CHROM", "POS")], by = "SNP") %>% filter(!is.na(p))
    w <- wza(d$CHROM, d$POS, d$p, d$af) %>% mutate(set = s, trait = tr)
    top <- w %>% filter(emp_p <= CFG$wza_top)
    h <- hits %>% filter(set == s, trait == tr)
    top$overlaps_EGWAS <- if (nrow(top)) unname(mapply(function(ch, st) any(h$CHROM == ch & h$POS >= st & h$POS < st + CFG$wza_window),
                                                     top$chr, top$start)) else logical(0)
    out[[paste(s, tr)]] <- top
  }
}
s9 <- bind_rows(out); write.csv(s9, tab("TableS9_WZA_top_windows.csv"), row.names = FALSE)
if (nrow(s9)) print(s9 %>% group_by(set) %>% summarise(top_windows = n(), distinct_windows = n_distinct(win),
                                         overlapping_EGWAS = sum(as.logical(overlaps_EGWAS)))) else msg("No WZA windows in the top fraction")
