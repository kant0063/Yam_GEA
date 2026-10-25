################################################################################
# 01_build_master_data.R                                   Analysis Plan: A0
# Builds ONE aligned master object used by every later script:
#   dosage matrix, SNP map, metadata (groups, species, admixture, country, coords),
#   environmental table (+ VIF-pruned set), genetic PCs and kinship.
# Replaces the per-script VCF reading / sample matching in all original scripts,
# and the PC step of Yam_PCA.Rmd (Fig. 1A can still be drawn from Yam_PCA.Rmd).
# Fixes: explicit GT->dosage; environmental rows matched to genotype rows by ID
# (original Yam_RDA.R used unmatched env_data rows); latitude/longitude excluded
# from environmental variables regardless of capitalisation (original
# Yam_egwas_rrblup_V2.R excluded only lower-case names).
# Run:  Rscript 01_build_master_data.R
#       YAM_GENO_SET=relaxed Rscript 01_build_master_data.R   (after 03)
################################################################################
source("00_config_and_helpers.R")

vcf <- switch(CFG$geno_set, strict = CFG$vcf_strict, relaxed = CFG$vcf_relaxed)
msg("Reading", vcf)
gd <- read_dosage(vcf); G <- gd$G; snp <- gd$snp
msg("Genotypes:", nrow(G), "samples x", ncol(G), "SNPs; missing =", sum(is.na(G)))
if (CFG$geno_set == "strict" && any(is.na(G))) warning("Strict set should have no missing data")

# ---- metadata --------------------------------------------------------------
meta <- read.csv(CFG$meta_csv, stringsAsFactors = FALSE)
stopifnot(all(c("SampleID", "Population.assignment") %in% names(meta)))
if (file.exists(CFG$accession_csv)) {
  acc <- read.csv(CFG$accession_csv, stringsAsFactors = FALSE)       # SampleID, Field_species, Country
  meta <- left_join(meta, acc, by = "SampleID")
} else msg("NOTE:", CFG$accession_csv, "not found - Field_species/Country left NA (needed for Table 1)")

pa <- meta$Population.assignment
meta$Hybrid <- grepl(" x |hybrid|admix", pa, ignore.case = TRUE)   # "Hybrid" label or "A x B" labels
meta$Species_group <- case_when(meta$Hybrid ~ "hybrid",
                                grepl("rotundata", pa) ~ "rotundata",
                                grepl("praehensilis", pa) ~ "praehensilis",
                                grepl("abyssinica", pa) ~ "abyssinica", TRUE ~ "unassigned")
if (any(meta$Species_group == "unassigned")) {
  msg("WARNING:", sum(meta$Species_group == "unassigned"), "samples have a Population.assignment that names no species;",
      "they stay in the full panel but are excluded from within-species analyses:")
  print(table(meta$Population.assignment[meta$Species_group == "unassigned"], useNA = "ifany"))
}
meta$Lineage <- ifelse(grepl("Cameroon", pa), "Cameroonian", ifelse(grepl("Western", pa), "Western", NA))

# admixture (TESS3 K = 7 Q matrix written by Yam_spatial_population_structure.R)
if (file.exists(CFG$qmatrix_csv)) {
  Q <- read.csv(CFG$qmatrix_csv, stringsAsFactors = FALSE)
  idcol <- intersect(c("SampleID", "Sample", "ID", "Accession"), names(Q))
  qnum  <- Q[, vapply(Q, is.numeric, logical(1)), drop = FALSE]
  qnum  <- qnum[, !names(qnum) %in% c("Longitude", "Latitude", "lon", "lat"), drop = FALSE]
  Qmax  <- apply(qnum, 1, max)
  ids   <- if (length(idcol)) Q[[idcol[1]]] else rownames(G)   # assumes VCF order if no ID column
  meta$Qmax <- Qmax[match(meta$SampleID, ids)]
} else { msg("NOTE: no Q matrix; NonAdmixed defined from Population.assignment only"); meta$Qmax <- NA }
meta$NonAdmixed <- !meta$Hybrid & meta$Species_group != "unassigned" & (is.na(meta$Qmax) | meta$Qmax >= CFG$q_pure)

# ---- environment -----------------------------------------------------------
env <- read.csv(CFG$env_csv, stringsAsFactors = FALSE, check.names = FALSE)
names(env)[1] <- "SampleID"
lat <- grep("^lat", names(env), ignore.case = TRUE, value = TRUE)[1]
lon <- grep("^lon", names(env), ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(lat), !is.na(lon))
env_vars <- setdiff(names(env)[vapply(env, is.numeric, logical(1))], c(lat, lon))
msg("Environmental variables (", length(env_vars), "):", paste(env_vars, collapse = ", "))

# ---- align everything to genotype order -------------------------------------
ids <- Reduce(intersect, list(rownames(G), meta$SampleID, env$SampleID))
msg("Samples in VCF, metadata and env table:", length(ids), "(VCF:", nrow(G), ")")
G <- G[ids, ]; meta <- meta[match(ids, meta$SampleID), ]; env <- env[match(ids, env$SampleID), ]
meta$Latitude <- env[[lat]]; meta$Longitude <- env[[lon]]
meta$SiteID <- as.integer(factor(paste(round(meta$Latitude, 4), round(meta$Longitude, 4))))
E <- as.matrix(env[, env_vars]); rownames(E) <- ids
if (any(is.na(E))) { msg("WARNING:", sum(is.na(E)), "missing env values - imputed with column means (report!)")
  E <- apply(E, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x }) }
E <- E[, apply(E, 2, sd) > 0, drop = FALSE]

# VIF pruning (A0) - all variables, and bioclim-only set used for genomic offset (A8)
vp_all <- vif_prune(scale(E))
bio <- grep("^bio", colnames(E), ignore.case = TRUE, value = TRUE)
vp_bio <- vif_prune(scale(E[, bio]))
msg("VIF <", CFG$vif_max, "kept", length(vp_all$keep), "of", ncol(E), ":", paste(vp_all$keep, collapse = ", "))

# ---- genetic PCs and kinship -----------------------------------------------
Gf <- maf_filter(G); keepT <- thin_snps(snp[match(colnames(Gf), snp$SNP), ])
pca <- prcomp(Gf[, keepT], center = TRUE, scale. = FALSE)
PCs <- pca$x[, 1:10]; pc_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
msg("PC1-3 variance %:", paste(round(pc_var[1:3], 2), collapse = ", "), " (manuscript reports 2.66, 2.43)")
K <- kinship(G)

M <- list(G = G, snp = snp, meta = meta, E = E, env_vars = colnames(E),
          vif_keep = vp_all$keep, vif_keep_bio = vp_bio$keep, PCs = PCs, pc_var = pc_var,
          K = K, thin = thin_snps(snp), cfg = CFG)
saveRDS(M, res("master.rds"))
write.csv(meta, tab("sample_metadata_aligned.csv"), row.names = FALSE)
write.csv(data.frame(variable = colnames(E), VIF_retained = colnames(E) %in% vp_all$keep,
                     bioclim_offset_set = colnames(E) %in% vp_bio$keep), tab("TableS1_environmental_variables_VIF.csv"), row.names = FALSE)
print(table(meta$Population.assignment, meta$NonAdmixed, dnn = c("group", "non-admixed")))
msg("Saved", res("master.rds"))
