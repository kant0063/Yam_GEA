################################################################################
# 06_Fst_outliers_revised.R                                Analysis Plan: A6c
# REPLACES the FST sections of integrated_population_genomics_analysis.R and
# Yam_Fst_outlier_analysis.R (ROH and private-allele sections of the integrated
# script can be kept as they are).
# Why: the originals call hierfstat::betas() on 0/1/2 dosages (mean-imputed, so
# not integers). betas() expects allele-coded genotypes (e.g. 11/12/22) and
# returns Weir & Hill betas, not Weir & Cockerham FST - check the old values.
# This script computes per-locus Weir & Cockerham (1984) FST directly and makes
# the contrast explicit: (1) global across the 7 genetic groups (what Fig. 2C
# showed), (2) pairwise among the four non-admixed groups.
# Outputs: Fig. 2C (revised), Fig. S4 (revised), outlier tables.
################################################################################
source("00_config_and_helpers.R"); M <- load_master(); meta <- M$meta; G <- M$G

grp <- meta$Population.assignment
glob <- wc_fst(G, grp)
fst <- data.frame(M$snp[match(colnames(G), M$snp$SNP), c("SNP", "CHROM", "POS")], FST_global = glob$fst)
msg("Genome-wide WC FST across", length(unique(grp)), "groups =", round(glob$global, 3))

# private alleles: alt allele present in exactly one group (strict set)
present <- sapply(split(seq_along(grp), grp), function(i) colSums(G[i, , drop = FALSE], na.rm = TRUE) > 0)
fst$private_in <- ifelse(rowSums(present) == 1, colnames(present)[max.col(present)], NA)

# pairwise contrasts among non-admixed groups
pure <- meta$NonAdmixed; pg <- unique(grp[pure]); pairs <- combn(pg, 2, simplify = FALSE)
pw <- sapply(pairs, function(pp) { i <- grp %in% pp & pure; wc_fst(G[i, ], grp[i])$fst })
colnames(pw) <- sapply(pairs, function(pp) paste(gsub("D\\. ", "", pp), collapse = " vs "))
fst <- cbind(fst, pw)
pw_global <- sapply(pairs, function(pp) { i <- grp %in% pp & pure; wc_fst(G[i, ], grp[i])$global })
write.csv(data.frame(contrast = colnames(pw), WC_FST = pw_global), tab("FST_pairwise_genomewide.csv"), row.names = FALSE)

q <- quantile(fst$FST_global, c(.95, .99, .995), na.rm = TRUE)
fst$outlier <- cut(fst$FST_global, c(-Inf, q, Inf), labels = c("none", "95th", "99th", "99.5th"))
write.csv(fst, res("fst_per_locus_WC.csv"), row.names = FALSE)
write.csv(fst %>% filter(outlier != "none"), tab("FST_outliers_global.csv"), row.names = FALSE)
write.csv(fst %>% filter(outlier != "none") %>% count(CHROM, outlier), tab("FigS4_outliers_by_chromosome.csv"), row.names = FALSE)
msg("Outliers with private alleles:", sum(fst$outlier != "none" & !is.na(fst$private_in), na.rm = TRUE),
    "of", sum(fst$outlier != "none", na.rm = TRUE))

# Fig. 2C revised: contrast stated in the title
fst <- fst %>% arrange(CHROM, POS) %>% mutate(idx = row_number())
p <- plot_fst_by_chromosome(fst)
ggsave(fig("Fig2C_FST_global_WC.png"), p, width = 13, height = 4.5, dpi = 300, bg = "white")   # = composite panel A
