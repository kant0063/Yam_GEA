################################################################################
# 00_config_and_helpers.R
# Yam GEA revision (TPG-2026-07-0232) - shared configuration and helper functions
# Sourced by every numbered script (01-15). Edit the CONFIG block only.
#
# Run order:  01 -> 02 -> 03 (relaxed VCF; then rerun 01 with YAM_GENO_SET=relaxed)
#             -> 04 -> ... -> 17 -> 18 -> 19   (see README_run_order.md)
# Each script writes to results/<geno_set>/, figures/<geno_set>/ and tables/<geno_set>/
################################################################################

# Remove objects left in the R session that mask base functions used here (e.g. an older
# script 19 defined t() and r(), which broke t() in script 01 when run in the same session).
for (.nm in c("t", "c", "r", "q", "df")) if (exists(.nm, envir = globalenv(), inherits = FALSE) &&
    is.function(get(.nm, envir = globalenv()))) { rm(list = .nm, envir = globalenv()); message("Removed masking function '", .nm, "' from the R session") }
rm(.nm)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(data.table)
})

# ============================== CONFIG ========================================
CFG <- list(
  # ---- input files (same names as the original scripts) ----
  vcf_strict    = "Yam_no_missing.vcf",                  # 63,142 SNPs, zero missing
  vcf_full      = "Calling_ALL_Rotundata_Allc05.vcf.gz", # 3,570,940 SNPs
  vcf_relaxed   = "Yam_relaxed_imputed.vcf.gz",          # created by 03
  meta_csv      = "Yam_meta_data.csv",                   # SampleID, Population.assignment
  env_csv       = "Yam_bioclim_data_extracted_V3.csv",   # SampleID, lat, lon, 19 bio + 17 soil
  accession_csv = "Yam_accession_info.csv",              # OPTIONAL: SampleID, Field_species, Country (from the original Table S6, accession list)
  qmatrix_csv   = "tess3_qmatrix_K7.csv",                # written by Yam_spatial_population_structure.R
  koppen_hist   = "koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif",
  koppen_fut    = "koppen_geiger_tif/2071_2099/ssp585/koppen_geiger_0p00833333.tif",
  clim_dir      = "climate_data",                        # geodata download cache
  gff           = Sys.getenv("YAM_GFF", "TDr96_annotation.gff3"),  # D. rotundata TDr96 annotation (16); give a full path if
                                                                    # it is not in the scripts folder or its parent
  proteins_fa   = "TDr96_proteins.fa",                   # OPTIONAL: for BLAST in 16
  tair_db       = "TAIR10_pep",                          # OPTIONAL: BLAST db of TAIR10 proteins (16)

  # ---- which genotype set to use: "strict" (default) or "relaxed" (sensitivity, A2) ----
  geno_set      = Sys.getenv("YAM_GENO_SET", "strict"),

  # ---- decisions (Analysis Plan section 2) ----
  q_pure        = 0.80,   # max ancestry >= q_pure -> non-admixed
  min_n_within  = 20,     # skip within-species sets smaller than this
  n_pc          = 6,      # genetic PCs used as ancestry covariates
  vif_max       = 10,     # VIF threshold for environmental predictors
  fdr           = 0.05,
  min_ve_ratio  = 0.05,   # REML lower bound on Ve/Vu (h2 <= ~0.95); avoids boundary fits that make GLS unstable
  int_transform = TRUE,   # rank-based inverse-normal transform of EGWAS responses (removes outlier leverage)   # FDR for EGWAS / LFMM (original V2 script used 0.20)
  maf           = 0.05,
  thin_bp       = 7000,   # SNP thinning distance for RDA permutations (~LD half-decay)
  perms         = 999,
  rda_cores     = as.integer(Sys.getenv("YAM_RDA_CORES", "4")),    # parallel permutation workers (each needs memory)
  rda_max_snps  = as.integer(Sys.getenv("YAM_RDA_MAX_SNPS", "0")), # 0 = no cap; else random subset of thinned SNPs
  wza_window    = 10000,  # bp
  wza_top       = 0.001,  # top 0.1% of windows
  relaxed_missing = 0.20,
  K_lfmm_full   = 7,
  K_lfmm_within = c(rotundata = 3, praehensilis = 2),   # check with snmf cross-entropy (09)
  ssp           = "585",
  period        = "2081-2100",
  gcms          = c("ACCESS-ESM1-5","BCC-CSM2-MR","CanESM5","CNRM-CM6-1","IPSL-CM6A-LR","MIROC6"),
  clim_res      = 10,     # arc-minutes
  map_buffer    = 3,      # degrees around sampled extent
  core_csv      = "Yam_core_collection.csv",  # original 50-accession core (column SampleID); else Table S14 from 13
  core_reps     = 50,     # random size-matched training sets compared with the core
  cv_k          = 10,
  cv_reps       = 5,
  n_spatial_blocks = 6,
  heavy_traits  = c("IC1","IC2","IC3","bio_12","bio_17","PHIHOX"),  # traits for BayesB/RKHS/RF
  bglr_iter     = 6000, bglr_burn = 1000,
  power_reps    = 100,
  power_pve     = c(0.01, 0.02, 0.05, 0.10, 0.20),
  gemma_bin     = Sys.getenv("YAM_GEMMA", "/Users/mikey/Downloads/GEMMA/bin/gemma"),  beagle_jar = "beagle.jar", bcftools = "bcftools",
  roh_min_snp   = 20, roh_max_opp = 1, roh_max_miss = 1,   # same as original integrated script
  roh_min_bp    = 1e6, roh_max_gap = 2.5e5,                 # ROH >= 1 Mb, gaps <= 250 kb
  rarefy_genes  = 10,     # gene copies for rarefied private-allele richness (<= 2 x smallest group)
  ld_max_dist   = 100000, ld_anchors = 3000,
  gene_window   = 10000,  # bp each side of a marker for candidate genes (or LD half-decay, if larger)
  seed          = 2026
)
set.seed(CFG$seed)
# look for the GFF in the scripts folder, then the parent folder
if (!file.exists(CFG$gff) && file.exists(file.path("..", CFG$gff))) CFG$gff <- file.path("..", CFG$gff)

# names used for key variables (matched case-insensitively, with or without "_")
VAR_PRECIP <- "bio_12"; VAR_PH <- "PHIHOX"

# ============================== OUTPUT PATHS ==================================
RES <- file.path("results", CFG$geno_set); FIG <- file.path("figures", CFG$geno_set)
TAB <- file.path("tables", CFG$geno_set)
for (d in c(RES, FIG, TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
res <- function(x) file.path(RES, x); fig <- function(x) file.path(FIG, x); tab <- function(x) file.path(TAB, x)
need <- function(pk) { miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Install packages: ", paste(miss, collapse = ", ")) }
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# ============================== GENOTYPES =====================================
# Explicit GT -> alt-allele dosage (0/1/2). vcfR::extract.gt(as.numeric = TRUE) is not a
# reliable GT->dosage conversion, so the revised scripts never use it for GT.
gt_to_dosage <- function(gt) {
  g <- gsub("\\|", "/", gt)
  d <- matrix(NA_real_, nrow(g), ncol(g), dimnames = dimnames(g))
  d[g == "0/0"] <- 0; d[g %in% c("0/1", "1/0")] <- 1; d[g == "1/1"] <- 2
  d
}
read_dosage <- function(vcf_file) {
  need("vcfR")
  v <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
  keep <- nchar(v@fix[, "REF"]) == 1 & nchar(v@fix[, "ALT"]) == 1   # biallelic SNPs only
  v <- v[keep, ]
  G <- base::t(gt_to_dosage(vcfR::extract.gt(v, element = "GT")))          # samples x SNPs
  snp <- data.frame(SNP = paste(v@fix[, "CHROM"], v@fix[, "POS"], sep = "_"),
                    CHROM = v@fix[, "CHROM"], POS = as.numeric(v@fix[, "POS"]),
                    REF = v@fix[, "REF"], ALT = v@fix[, "ALT"], stringsAsFactors = FALSE)
  colnames(G) <- snp$SNP
  list(G = G, snp = snp)
}
load_master <- function() {
  f <- Sys.getenv("YAM_MASTER", res("master.rds"))   # 05.5 points this at a light, thinned copy if (!file.exists(f)) stop("Run 01_build_master_data.R first (geno_set = ", CFG$geno_set, ")")
  readRDS(f)
}
# sample subsets used throughout: full panel and within-species (non-admixed only)
# %in% (not ==) so that samples with NA species group / admixture are dropped, never returned as NA IDs
subsets <- function(M) {
  ok <- M$meta$NonAdmixed %in% TRUE
  list(full         = M$meta$SampleID[!is.na(M$meta$SampleID)],
       rotundata    = M$meta$SampleID[M$meta$Species_group %in% "rotundata"    & ok],
       praehensilis = M$meta$SampleID[M$meta$Species_group %in% "praehensilis" & ok])
}
# world outline layer: polygons from maps::map_data("world") drawn with geom_path.
# Works with every ggplot2 version (borders()/annotation_borders() changed in ggplot2 4.0 and
# failed inside faceted maps, which produced a blank Fig. S17).
world_borders <- function(xlim = NULL, ylim = NULL, colour = "grey50", linewidth = 0.3, ...) {
  need("maps"); w <- ggplot2::map_data("world")
  if (!is.null(xlim)) w <- w[w$group %in% unique(w$group[w$long >= xlim[1] - 5 & w$long <= xlim[2] + 5]), ]
  if (!is.null(ylim)) w <- w[w$group %in% unique(w$group[w$lat >= ylim[1] - 5 & w$lat <= ylim[2] + 5]), ]
  ggplot2::geom_path(data = w, ggplot2::aes(x = long, y = lat, group = group), inherit.aes = FALSE,
                     colour = colour, linewidth = linewidth)
}
maf_filter <- function(G, maf = CFG$maf) { f <- colMeans(G, na.rm = TRUE) / 2; G[, pmin(f, 1 - f) >= maf, drop = FALSE] }
thin_snps <- function(snp, bp = CFG$thin_bp) {     # keep first SNP, then next >= bp away (per chromosome)
  keep <- logical(nrow(snp)); o <- order(snp$CHROM, snp$POS); last_chr <- ""; last_pos <- -Inf
  for (i in o) { if (snp$CHROM[i] != last_chr || snp$POS[i] - last_pos >= bp) {
      keep[i] <- TRUE; last_chr <- snp$CHROM[i]; last_pos <- snp$POS[i] } }
  keep
}
kinship <- function(G) { need("rrBLUP"); rrBLUP::A.mat(G - 1) }   # VanRaden-type, input coded -1/0/1

# ============================== ENVIRONMENT ===================================
vif_prune <- function(E, thr = CFG$vif_max) {          # iterative removal of the highest-VIF variable
  E <- as.data.frame(E); repeat {
    v <- diag(solve(cor(E))); if (max(v) < thr || ncol(E) <= 2) break
    E <- E[, -which.max(v), drop = FALSE] }
  list(keep = colnames(E), vif = diag(solve(cor(E))))
}
find_var <- function(nms, target) {   # matches bio_12 / bio12 / BIO12 ...
  key <- function(x) gsub("[^a-z0-9]", "", tolower(x)); i <- which(key(nms) == key(target)); if (length(i)) nms[i[1]] else NA
}

# ============================== STATISTICS ====================================
# Weir & Cockerham (1984) FST per biallelic locus from dosage matrix G (samples x SNPs)
wc_fst <- function(G, pop) {
  pop <- droplevels(as.factor(pop)); L <- levels(pop); r <- length(L)
  ni <- sapply(L, function(l) colSums(!is.na(G[pop == l, , drop = FALSE])))
  pi <- sapply(L, function(l) colMeans(G[pop == l, , drop = FALSE], na.rm = TRUE) / 2)
  hi <- sapply(L, function(l) colMeans(G[pop == l, , drop = FALSE] == 1, na.rm = TRUE))
  if (is.null(dim(ni))) { ni <- base::t(ni); pi <- base::t(pi); hi <- base::t(hi) }
  nbar <- rowSums(ni) / r
  nc   <- (r * nbar - rowSums(ni^2) / (r * nbar)) / (r - 1)
  pbar <- rowSums(ni * pi) / (r * nbar)
  s2   <- rowSums(ni * (pi - pbar)^2) / ((r - 1) * nbar)
  hbar <- rowSums(ni * hi) / (r * nbar)
  a <- nbar / nc * (s2 - 1 / (nbar - 1) * (pbar * (1 - pbar) - (r - 1) / r * s2 - hbar / 4))
  b <- nbar / (nbar - 1) * (pbar * (1 - pbar) - (r - 1) / r * s2 - (2 * nbar - 1) / (4 * nbar) * hbar)
  c <- hbar / 2
  list(fst = a / (a + b + c), global = sum(a, na.rm = TRUE) / sum(a + b + c, na.rm = TRUE))
}

# EMMAX-style single-SNP scan (P3D): one REML fit per trait, then GLS for every SNP.
# Returns beta (per alt allele), SE, p-value, allele frequency and PVE. Covariates X optional.
# Design matrix with intercept; drops constant, NA-only and collinear covariates so that
# rrBLUP::mixed.solve never sees a rank-deficient X ("X not full rank").
full_rank_X <- function(X, n) {
  X <- if (is.null(X)) matrix(1, n, 1) else cbind(1, as.matrix(X))
  X[is.na(X)] <- 0
  q <- qr(X); keep <- sort(q$pivot[seq_len(q$rank)])
  if (!1 %in% keep) keep <- c(1, keep[-length(keep)])
  if (length(keep) < ncol(X)) msg("Dropped", ncol(X) - length(keep), "constant/collinear covariate(s)")
  X[, keep, drop = FALSE]
}
emmax_scan <- function(y, G, K, X = NULL, fit = NULL) {
  need("rrBLUP"); n <- length(y)
  X <- full_rank_X(X, n)
  if (is.null(fit)) fit <- rrBLUP::mixed.solve(y, X = X, K = K, bounds = c(CFG$min_ve_ratio %||% 0.05, 1e9))
  G <- as.matrix(G); storage.mode(G) <- "double"
  V <- fit$Vu * K + fit$Ve * diag(n)
  e <- eigen(V, symmetric = TRUE); Vih <- e$vectors %*% (base::t(e$vectors) / sqrt(e$values))
  Q <- qr(Vih %*% X)
  yr <- qr.resid(Q, Vih %*% y); Gr <- qr.resid(Q, Vih %*% G)
  ssg <- colSums(Gr^2); beta <- colSums(Gr * as.vector(yr)) / ssg
  df <- n - ncol(X) - 1; rss <- sum(yr^2) - beta^2 * ssg
  se <- sqrt(pmax(rss, 0) / df / ssg); p <- 2 * pt(-abs(beta / se), df)
  af <- colMeans(G, na.rm = TRUE) / 2
  # PVE = share of residual (post-covariate, GLS) variance explained by the SNP: t^2/(t^2 + df); bounded in [0, 1]
  tt <- beta / se; pve <- tt^2 / (tt^2 + df)
  bad <- ssg < 1e-10; beta[bad] <- se[bad] <- p[bad] <- pve[bad] <- NA
  data.frame(SNP = colnames(G), beta = beta, se = se, p = p, af = af, PVE = pve,
             h2_null = fit$Vu / (fit$Vu + fit$Ve), row.names = NULL)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
# rank-based inverse-normal transform (Blom offset)
rint <- function(x) { r <- rank(x, na.last = "keep", ties.method = "average"); qnorm((r - 0.375) / (sum(!is.na(x)) + 0.25)) }
lambda_gc <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }

# Weighted-Z analysis (Booker et al. 2024): per-SNP one-sided Z from GEA p-values, weighted by pq,
# combined within windows; Z_W standardised within bins of SNP count to remove window-size effects.
wza <- function(chr, pos, p, af, window = CFG$wza_window) {
  d <- data.frame(win = sprintf("%s:%.0f", chr, floor(pos / window) * window), chr = chr,
                  start = floor(pos / window) * window, z = qnorm(pmin(pmax(1 - p, 1e-15), 1 - 1e-15)),
                  w = af * (1 - af))
  W <- d %>% group_by(win, chr, start) %>%
    summarise(nSNP = n(), ZW = sum(w * z) / sqrt(sum(w^2)), .groups = "drop") %>%
    filter(is.finite(ZW)) %>%
    mutate(bin = cut(nSNP, unique(quantile(nSNP, 0:10 / 10)), include.lowest = TRUE)) %>%
    group_by(bin) %>% mutate(ZW_std = (ZW - mean(ZW)) / sd(ZW)) %>% ungroup() %>%
    mutate(emp_p = rank(-ZW_std) / n())
  W
}

# Multiple matrix regression with randomization (Wang 2013)
mmrr <- function(Y, Xlist, nperm = CFG$perms) {
  lt <- lower.tri(Y); y <- as.numeric(scale(Y[lt]))
  X <- sapply(Xlist, function(m) as.numeric(scale(m[lt]))); colnames(X) <- names(Xlist)
  s0 <- summary(lm(y ~ X)); t0 <- s0$coefficients[-1, 3]; F0 <- s0$fstatistic[1]
  tp <- rep(1, length(t0)); Fp <- 1; n <- nrow(Y)
  for (i in seq_len(nperm)) { o <- sample(n); s <- summary(lm(as.numeric(scale(Y[o, o][lt])) ~ X))
    tp <- tp + (abs(s$coefficients[-1, 3]) >= abs(t0)); Fp <- Fp + (s$fstatistic[1] >= F0) }
  data.frame(term = names(Xlist), coef = s0$coefficients[-1, 1], t = t0, p = tp / (nperm + 1),
             R2 = s0$r.squared, F = F0, F_p = Fp / (nperm + 1), row.names = NULL)
}

# ExDet (Mesgaran et al. 2014): NT1 < 0 = univariate novelty; NT2 > 1 = combinatorial novelty
exdet <- function(ref, new) {
  ref <- as.matrix(ref); new <- as.matrix(new)
  mn <- apply(ref, 2, min); mx <- apply(ref, 2, max); rg <- mx - mn; rg[rg == 0] <- 1
  nt1 <- rowSums(pmin(sweep(new, 2, mn, "-") / rep(rg, each = nrow(new)),
                      sweep(-new, 2, -mx, "-") / rep(rg, each = nrow(new)), 0))
  mu <- colMeans(ref); S <- cov(ref)
  dref <- mahalanobis(ref, mu, S); dnew <- mahalanobis(new, mu, S)
  nt2 <- dnew / max(dref)
  data.frame(NT1 = nt1, NT2 = nt2, analog = nt1 == 0 & nt2 <= 1)
}


# ============================== FIGURE HELPERS ================================
# Fig. S10: within-group GEAV heatmap. Rows labelled with accession names (candidate parents
# marked with *), ordered within each group by the first within-group genetic PC; a side bar
# shows k-means sub-clusters (k = 2) on within-group PCs, so residual structure inside a group
# (e.g. two D. abyssinica sub-groups) is visible and quantified (Table S16b).
plot_geav_heatmap <- function(out, M, cand_ids = character(0), cap = 3) {
  meta <- M$meta
  namecol <- intersect(c("Name", "Name.x", "Accession", "Accession_name"), names(meta))
  lab <- if (length(namecol)) as.character(meta[[namecol[1]]][match(out$SampleID, meta$SampleID)]) else out$SampleID
  lab[is.na(lab) | lab == ""] <- out$SampleID[is.na(lab) | lab == ""]
  lab <- ifelse(out$SampleID %in% cand_ids, paste0(lab, " *"), lab)
  # within-group genetic PC1 and k-means sub-clusters
  sub <- data.frame(SampleID = out$SampleID, group = out$group, pc1 = NA_real_, cluster = NA_character_, sep_ratio = NA_real_)
  for (g in unique(out$group)) {
    i <- which(out$group == g); if (length(i) < 8) next
    Gs <- maf_filter(M$G[out$SampleID[i], M$thin, drop = FALSE]); Gs[is.na(Gs)] <- 0
    pc <- prcomp(Gs)$x[, 1:2]; set.seed(CFG$seed); km <- kmeans(pc, 2, nstart = 25)
    sub$pc1[i] <- pc[, 1]; sub$cluster[i] <- paste0(g, " - ", c("a", "b")[km$cluster])
    sub$sep_ratio[i] <- km$betweenss / km$tot.withinss          # >1: clusters more distinct than diffuse
  }
  hm <- out %>% tidyr::pivot_longer(-c(SampleID, group), names_to = "variable", values_to = "z") %>%
    dplyr::group_by(variable) %>% dplyr::mutate(z = as.numeric(scale(z))) %>% dplyr::ungroup()
  vt <- function(v) { k <- suppressWarnings(as.integer(gsub("\\D", "", v)))
    ifelse(grepl("^IC", v), "IC", ifelse(grepl("^bio", v, ignore.case = TRUE), ifelse(!is.na(k) & k >= 12, "precipitation", "temperature"), "soil")) }
  vv <- unique(hm$variable); vord <- vv[order(factor(vt(vv), c("temperature", "precipitation", "soil", "IC")), vv)]
  ord <- order(out$group, sub$pc1)
  hm$variable <- factor(hm$variable, levels = c("sub-cluster", vord))
  hm$row <- factor(hm$SampleID, levels = out$SampleID[ord])
  bar <- data.frame(SampleID = out$SampleID, group = out$group, variable = factor("sub-cluster", levels = levels(hm$variable)),
                    row = factor(out$SampleID, levels = out$SampleID[ord]), cl = substr(sub$cluster, nchar(sub$cluster), nchar(sub$cluster)))
  p <- ggplot(hm, aes(variable, row)) + geom_tile(aes(fill = z)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-cap, cap), oob = scales::squish,
                         name = paste0("Within-group\nGEAV (z,\ncapped at \u00b1", cap, ")")) +
    ggnewscale_safe(bar) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = setNames(lab, out$SampleID)) +
    labs(x = NULL, y = NULL) + theme_minimal() +
    theme(axis.text.y = element_text(size = 3.2), axis.text.x = element_text(angle = 90, size = 6, hjust = 1, vjust = .5),
          strip.text.y = element_text(angle = 0, size = 7), panel.grid = element_blank())
  attr(p, "subclusters") <- sub
  p
}
# sub-cluster side bar drawn as outlined tiles (no second fill scale needed)
ggnewscale_safe <- function(bar) {
  list(geom_tile(data = bar[!is.na(bar$cl) & bar$cl == "a", ], aes(variable, row), fill = "grey35", inherit.aes = FALSE),
       geom_tile(data = bar[!is.na(bar$cl) & bar$cl == "b", ], aes(variable, row), fill = "grey80", inherit.aes = FALSE))
}


# Fig. 2A: per-SNP global FST along the genome, x-axis labelled by pseudomolecule (cumulative bp),
# alternating grey by pseudomolecule, private-allele SNPs as green triangles, percentile lines.
chrom_number <- function(x) {          # "BDMI01000021.1" -> 21; "chr5" -> 5
  n <- suppressWarnings(as.numeric(sub("^.*?(\\d+)$", "\\1", sub("\\.\\d+$", "", x), perl = TRUE)))
  ifelse(!is.na(n) & n >= 1000, n %% 1000, n)
}
plot_fst_by_chromosome <- function(fs, point_size = .35) {
  fs <- fs[!is.na(fs$FST_global), ]
  chr_num <- chrom_number(fs$CHROM)
  ord <- unique(fs$CHROM[order(chr_num, fs$CHROM)])
  fs$CHROM <- factor(fs$CHROM, levels = ord)
  len <- tapply(fs$POS, fs$CHROM, max); off <- c(0, cumsum(as.numeric(len))[-length(len)]); names(off) <- names(len)
  fs$x <- fs$POS + off[as.character(fs$CHROM)]
  mid <- off + len / 2; lab <- chrom_number(names(len))
  qq <- quantile(fs$FST_global, c(.95, .99, .995), na.rm = TRUE)
  ggplot(fs, aes(x / 1e6, FST_global)) +
    geom_point(aes(colour = CHROM), size = point_size, show.legend = FALSE) +
    scale_colour_manual(values = rep(c("grey30", "grey70"), length.out = nlevels(fs$CHROM))) +
    geom_point(data = fs[!is.na(fs$private_in), ], colour = "darkgreen", shape = 17, size = point_size * 1.6) +
    geom_hline(yintercept = qq, linetype = 2, colour = c("gold2", "orange", "red")) +
    scale_x_continuous(breaks = mid / 1e6, labels = ifelse(is.na(lab), names(len), lab), expand = expansion(mult = .005)) +
    labs(x = "Pseudomolecule", y = expression(F[ST])) + theme_bw() +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 7))
}
# Fig. 2 composite: A FST by pseudomolecule; B rarefied and C raw private alleles; D FROH
fig2_composite <- function(fs, pa, fr) {
  need("patchwork")
  lev <- pa$group[order(pa$private_rarefied)]
  pa$group <- factor(pa$group, levels = lev); fr$group <- factor(fr$group, levels = lev)
  b <- ggplot(pa, aes(group, private_rarefied)) + geom_col(fill = "grey35") + coord_flip() +
    scale_y_continuous(n.breaks = 4) + labs(x = NULL, y = sprintf("Private alleles (rarefied, %d copies)", CFG$rarefy_genes)) + theme_bw()
  c <- ggplot(pa, aes(group, private_raw)) + geom_col(fill = "grey35") + coord_flip() +
    scale_y_continuous(n.breaks = 4) + labs(x = NULL, y = "Private alleles (raw counts)") + theme_bw() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  d <- ggplot(fr, aes(group, FROH)) + geom_boxplot(outlier.size = .8, width = .6) +
    stat_summary(fun = mean, geom = "point", shape = 23, fill = "white", size = 2) + coord_flip() +
    labs(x = NULL, y = expression(F[ROH])) + theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  patchwork::wrap_elements(full = plot_fst_by_chromosome(fs)) / (b | c | d) + patchwork::plot_layout(heights = c(1, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
}


# ---- panel store: scripts save the ggplot objects of multi-part figures; script 17 assembles them ----
save_panel <- function(p, name) {
  d <- res("panels"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  saveRDS(p, file.path(d, paste0(name, ".rds"))); invisible(p)
}
load_panel <- function(name, script = "") {
  f <- file.path(res("panels"), paste0(name, ".rds"))
  if (!file.exists(f)) stop("panel ", name, " not found - rerun ", script)
  readRDS(f)
}


# Fig. S6 panels from saved EGWAS p-values: Manhattan by pseudomolecule + QQ (same style as script 07)
egwas_manhattan_qq <- function(Ptr, snp) {
  d <- Ptr[is.finite(Ptr$p) & Ptr$p > 0, ]
  d <- merge(d, snp[, c("SNP", "CHROM", "POS")], by = "SNP")
  d$chr <- chrom_number(d$CHROM); d <- d[order(d$chr, d$POS), ]; d$idx <- seq_len(nrow(d))
  d$q <- p.adjust(d$p, "BH"); sig <- d$q < CFG$fdr; thr <- if (any(sig)) max(d$p[sig]) else NA
  ax <- aggregate(idx ~ chr, d, function(v) (min(v) + max(v)) / 2)
  pm <- ggplot(d, aes(idx, -log10(p), colour = factor(chr %% 2))) + geom_point(size = .4, show.legend = FALSE) +
    scale_colour_manual(values = c("grey60", "black")) +
    scale_x_continuous(breaks = ax$idx, labels = ax$chr, expand = expansion(mult = .01)) +
    { if (!is.na(thr)) geom_hline(yintercept = -log10(thr), colour = "red", linetype = 2) } +
    labs(x = "Pseudomolecule", y = expression(-log[10](p))) + theme_bw() +
    theme(axis.text.x = element_text(size = 6), panel.grid.minor = element_blank())
  qq <- data.frame(e = -log10(ppoints(nrow(d))), o = -log10(sort(d$p)))
  pq <- ggplot(qq, aes(e, o)) + geom_point(size = .5) + geom_abline(slope = 1, intercept = 0, colour = "red") +
    labs(x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p))) + theme_bw()
  list(manhattan = pm, qq = pq, lambda = lambda_gc(d$p), n_sig = sum(sig))
}
