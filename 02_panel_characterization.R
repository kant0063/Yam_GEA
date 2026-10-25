################################################################################
# 02_panel_characterization.R                              Analysis Plan: A1
# Table 1 (composition, sites, hybrids, He, BIO12/pH ranges), Fig. S1 (SNPs per
# chromosome + MAF spectrum per group), Table S2/S3 part 1 (effect of zero-missing
# filter per group, from the full 3.57 M SNP VCF).   New script (no original).
################################################################################
source("00_config_and_helpers.R"); M <- load_master(); G <- M$G; meta <- M$meta

# ---- Table 1 ----------------------------------------------------------------
v_pr <- find_var(colnames(M$E), VAR_PRECIP); v_ph <- find_var(colnames(M$E), VAR_PH)
grp_he <- function(ids) { p <- colMeans(G[ids, , drop = FALSE], na.rm = TRUE) / 2; n <- 2 * length(ids)
  mean(2 * p * (1 - p) * n / (n - 1)) }
sp_col <- if ("Field_species" %in% names(meta) && any(!is.na(meta$Field_species))) "Field_species" else "Species_group"
t1 <- meta %>% mutate(Species = .data[[sp_col]]) %>% group_by(Species) %>%
  summarise(n = n(), Unique_sites = n_distinct(SiteID), Hybrids_K7 = sum(Hybrid),
            BIO12_range = paste(round(range(M$E[SampleID, v_pr])), collapse = "-"),
            pH_range = paste(round(range(M$E[SampleID, v_ph]), 2), collapse = "-"), .groups = "drop")
if ("Country" %in% names(meta)) t1 <- left_join(t1, meta %>% mutate(Species = .data[[sp_col]]) %>%
     count(Species, Country) %>% pivot_wider(names_from = Country, values_from = n, values_fill = 0), by = "Species")
he <- meta %>% group_by(Population.assignment) %>% summarise(n = n(), He = grp_he(SampleID),
       BIO12_range = paste(round(range(M$E[SampleID, v_pr])), collapse = "-"),
       pH_range = paste(round(range(M$E[SampleID, v_ph]), 2), collapse = "-"), .groups = "drop")
write.csv(t1, tab("Table1_panel_composition.csv"), row.names = FALSE)
write.csv(he, tab("Table1_groups_He_ranges.csv"), row.names = FALSE)
msg("Unique sites:", n_distinct(meta$SiteID), "; hybrids:", sum(meta$Hybrid)); print(he)

# ---- Fig. S1: SNPs per chromosome + folded MAF per group --------------------
chr_n <- count(M$snp, CHROM)
maf_g <- bind_rows(lapply(split(meta$SampleID, meta$Population.assignment), function(ids) {
  p <- colMeans(G[ids, , drop = FALSE], na.rm = TRUE) / 2; data.frame(maf = pmin(p, 1 - p)) }), .id = "group")
p1 <- ggplot(chr_n, aes(CHROM, n)) + geom_col(fill = "grey40") + coord_flip() + labs(x = NULL, y = "SNPs") + theme_bw()
p2 <- ggplot(maf_g, aes(maf)) + geom_histogram(bins = 25) + facet_wrap(~group, scales = "free_y") +
  labs(x = "Minor allele frequency (within group)", y = "SNPs") + theme_bw()
save_panel(p1, "FigS1A"); save_panel(p2, "FigS1B")   # composite Figure S1 is built by script 17
ggsave(fig("FigS1a_SNPs_per_chromosome.png"), p1, width = 5, height = 6, dpi = 300)
ggsave(fig("FigS1b_MAF_by_group.png"), p2, width = 10, height = 7, dpi = 300)
write.csv(chr_n, tab("FigS1_snps_per_chromosome.csv"), row.names = FALSE)

# ---- Table S2/S3 part 1: which groups drive the zero-missing filter -------------
# Uses SNPRelate on the full VCF (3.57 M SNPs): per-group missing rate per SNP.
if (CFG$geno_set == "strict" && file.exists(CFG$vcf_full) && requireNamespace("SNPRelate", quietly = TRUE)) {
  gds <- res("full_calls.gds")
  if (!file.exists(gds)) SNPRelate::snpgdsVCF2GDS(CFG$vcf_full, gds, method = "biallelic.only")
  g <- SNPRelate::snpgdsOpen(gds); on.exit(SNPRelate::snpgdsClose(g))
  s_ids <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "sample.id"))
  tot <- length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "snp.id")))
  s7 <- bind_rows(lapply(split(meta$SampleID, meta$Population.assignment), function(ids) {
    ids <- intersect(ids, s_ids)
    mr <- SNPRelate::snpgdsSNPRateFreq(g, sample.id = ids, with.id = TRUE)$MissingRate
    data.frame(n = length(ids), pct_SNPs_with_missing = 100 * mean(mr > 0), mean_missing_rate = mean(mr)) }), .id = "group")
  s7$total_SNPs <- tot; s7$retained_strict <- ncol(G)
  write.csv(s7, tab("TableS2_missing_filter_by_group.csv"), row.names = FALSE); print(s7)
} else msg("Skipping Table S2 (needs strict set, full VCF and SNPRelate)")
