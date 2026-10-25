################################################################################
# 14_ROH_private_alleles.R                                 Analysis Plan: A1 / Fig. 2A-B
# REPLACES the ROH and private-allele sections of integrated_population_genomics_analysis.R.
# Problems in the original:
#   - consecutiveRUNS.run(..., ROHet = TRUE) detects runs of HETEROzygosity, so the
#     values plotted as FROH (Fig. 2B) were not runs of homozygosity -> ROHet = FALSE.
#   - private-allele counts were raw counts; they scale with sample size (groups range
#     from a few hybrids to 86 D. rotundata) -> add rarefied private allelic richness
#     (Kalinowski 2004) at CFG$rarefy_genes gene copies.
#   - v2: runs must be >= 1 Mb (CFG$roh_min_bp) with gaps <= 250 kb (CFG$roh_max_gap); with ~1 SNP per
#     7 kb, 20-SNP runs without a length floor inflated FROH (~0.5).
#   - PLINK ped/map no longer needed as input; written here from the master object.
# Outputs: Fig. 2A, Fig. 2B, froh_by_individual.csv, private_alleles_by_group.csv
################################################################################
source("00_config_and_helpers.R"); need("detectRUNS"); M <- load_master(); meta <- M$meta; G <- M$G
grp <- meta$Population.assignment

# ---- ROH (runs of HOMOzygosity) ---------------------------------------------
ped <- res("roh.ped"); map <- res("roh.map")
snp <- M$snp[match(colnames(G), M$snp$SNP), ]
chr_num <- as.integer(factor(snp$CHROM, levels = unique(snp$CHROM)))
fwrite(data.table(chr_num, snp$SNP, 0, snp$POS), map, sep = "\t", col.names = FALSE)
al <- function(d, ref, alt) { a1 <- ifelse(d == 2, alt, ref); a2 <- ifelse(d == 0, ref, alt)
  a1[is.na(d)] <- "0"; a2[is.na(d)] <- "0"; c(rbind(a1, a2)) }
pedrows <- base::t(sapply(seq_len(nrow(G)), function(i) al(G[i, ], snp$REF, snp$ALT)))
fwrite(data.table(grp_id = gsub(" ", "_", grp), iid = rownames(G), 0, 0, 0, -9, pedrows), ped, sep = " ", col.names = FALSE)
runs <- detectRUNS::consecutiveRUNS.run(genotypeFile = ped, mapFile = map, minSNP = CFG$roh_min_snp,
                                        ROHet = FALSE, maxOppRun = CFG$roh_max_opp, maxMissRun = CFG$roh_max_miss,
                                        minLengthBps = CFG$roh_min_bp %||% 1e6, maxGap = CFG$roh_max_gap %||% 2.5e5)
genome_len <- sum(tapply(snp$POS, snp$CHROM, max) - tapply(snp$POS, snp$CHROM, min))
froh <- data.frame(SampleID = rownames(G), group = grp) %>%
  left_join(runs %>% group_by(id) %>% summarise(n_roh = n(), roh_bp = sum(lengthBps)), by = c("SampleID" = "id")) %>%
  mutate(n_roh = coalesce(n_roh, 0L), roh_bp = coalesce(roh_bp, 0), FROH = roh_bp / genome_len)
write.csv(froh, tab("Fig2B_froh_by_individual.csv"), row.names = FALSE)
print(froh %>% group_by(group) %>% summarise(n = n(), mean_FROH = mean(FROH), max_FROH = max(FROH)))
p2b <- ggplot(froh, aes(group, FROH)) + geom_boxplot(outlier.size = .8) +
  stat_summary(fun = mean, geom = "point", shape = 23, fill = "white") + coord_flip() +
  labs(x = NULL, y = expression(F[ROH])) + theme_bw()
ggsave(fig("Fig2B_FROH.png"), p2b, width = 6, height = 4, dpi = 300)

# ---- private alleles: raw and rarefied ---------------------------------------
# allele counts per group: N_ij = copies of allele i in group j, N_j = genes sampled at the locus
cnt <- function(ids) { g <- G[ids, , drop = FALSE]; n <- 2 * colSums(!is.na(g)); alt <- colSums(g, na.rm = TRUE)
  list(n = n, alt = alt, ref = n - alt) }
C <- lapply(split(rownames(G), grp), cnt); L <- names(C)
raw <- sapply(L, function(j) { others <- setdiff(L, j)
  oalt <- Reduce(`+`, lapply(C[others], `[[`, "alt")); oref <- Reduce(`+`, lapply(C[others], `[[`, "ref"))
  sum(C[[j]]$alt > 0 & oalt == 0) + sum(C[[j]]$ref > 0 & oref == 0) })
g <- CFG$rarefy_genes
Qp <- function(Nij, Nj) ifelse(Nj >= g, 1 - exp(lchoose(Nj - Nij, g) - lchoose(Nj, g)), NA)   # P(allele seen in g genes)
rare <- sapply(L, function(j) { tot <- 0
  for (a in c("alt", "ref")) {
    qj <- Qp(C[[j]][[a]], C[[j]]$n)
    qo <- Reduce(`*`, lapply(setdiff(L, j), function(k) 1 - Qp(C[[k]][[a]], C[[k]]$n)))
    tot <- tot + sum(qj * qo, na.rm = TRUE) }
  tot })
pa <- data.frame(group = L, n = sapply(L, function(j) sum(grp == j)), private_raw = raw, private_rarefied = rare,
                 rarefied_to_genes = g)
write.csv(pa, tab("Fig2A_private_alleles_by_group.csv"), row.names = FALSE); print(pa)
# NOTE: the original Fig. 2A used the full 3.57 M SNP VCF; to repeat on it, run this block with
# YAM_GENO_SET=relaxed (or load the full VCF with read_dosage(CFG$vcf_full) on a large-memory node).
pa_long <- tidyr::pivot_longer(pa, c(private_rarefied, private_raw))
pa_long$name <- factor(pa_long$name, levels = c("private_rarefied", "private_raw"),
                       labels = c(sprintf("Rarefied (%d copies)", CFG$rarefy_genes), "Raw"))
p2a <- ggplot(pa_long, aes(reorder(group, value), value)) +
  geom_col() + facet_wrap(~name, scales = "free_x") + coord_flip() + labs(x = NULL, y = "Private alleles") + theme_bw()
ggsave(fig("Fig2A_private_alleles.png"), p2a, bg = "white", width = 8, height = 4, dpi = 300)
