################################################################################
# 18_assemble_tables_and_placeholders.R                    Analysis Plan: A10
# Builds Supplementary_Tables_S1-S19.xlsx (manuscript order) and and writes
# placeholder_values.md: the numbers needed to replace the highlighted [XX]
# placeholders in the revised manuscript and response letter.   New script.
################################################################################
source("00_config_and_helpers.R"); need("openxlsx")
# Read every table safely: empty files (e.g. a table with no significant results, written as
# a blank line) are listed in the index sheet instead of stopping the script.
safe_read <- function(f) {
  if (file.info(f)$size == 0) return(NULL)
  x <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || ncol(x) == 0 || (ncol(x) == 1 && all(names(x) %in% c("", "x", "V1")) && nrow(x) == 0)) NULL else x
}
# ==== Supplementary_Tables_S1-S19.xlsx: one sheet per supplementary table, in manuscript order ====
# (numbering of the ordered revision, 25 Oct 2026). Accessions are shown by name with the sequencing
# file name kept in the last column; Table S1 merges the definitions in env_variable_definitions.csv.
meta_all <- load_master()$meta
namecol <- intersect(c("Name", "Name.x", "Accession", "Accession_name"), names(meta_all))
acc_name <- function(ids) { n <- if (length(namecol)) as.character(meta_all[[namecol[1]]][match(ids, meta_all$SampleID)]) else ids
  ifelse(is.na(n) | n == "", ids, n) }
with_names <- function(x) {
  if (is.null(x) || !"SampleID" %in% names(x)) return(x)
  x <- data.frame(Accession = acc_name(x$SampleID), x[, setdiff(names(x), "SampleID"), drop = FALSE],
                  `Sequencing file` = x$SampleID, check.names = FALSE); x }
rdt <- function(f) { p <- if (file.exists(f)) f else tab(f); if (file.exists(p)) safe_read(p) else NULL }
captions <- c(
  "Table S1" = "Environmental variables: definitions, units and sources, and whether each was retained after VIF pruning (16 variables) and used in the genomic offset analysis (seven bioclimatic variables).",
  "Table S2" = "Effect of the zero-missing-data filter by genetic group: number of accessions, percentage of SNPs with at least one missing call, mean missing rate, total SNPs and SNPs retained by the zero-missing filter.",
  "Table S3" = "Comparison of key results between the zero-missing-data set (63,142 SNPs) and the relaxed dataset (≤ 20% missing, MAF ≥ 0.05, imputed with Beagle; 1,573,596 SNPs): Procrustes correlation of genetic PCs 1–4, the environment-only pRDA fraction, the number of full-panel EGWAS associations and the number of associated 10-kb windows shared across the full panel and both within-species sets.",
  "Table S4" = "Isolation by distance and by environment: MMRR coefficients and partial Mantel tests for the full panel and within species.",
  "Table S5" = "Description of Independent Components (ICs)",
  "Table S6" = "WZA windows in the top 0.1% within D. rotundata and D. praehensilis, with overlap with EGWAS candidates.",
  "Table S7" = "Genomic inflation factor (λ), null-model heritability, number of SNPs tested and number of significant associations (FDR 0.05) for every environmental variable and IC, for the full panel and within D. rotundata and D. praehensilis; “boundary” marks fits in which the variance component reached its bound.",
  "Table S8" = "LFMM significant SNPs (FDR 0.05) by environmental variable for the full panel (K = 7) and within D. rotundata (K = 3) and D. praehensilis (K = 2).",
  "Table S9" = "BSLMM posterior estimates of PVE, PGE and number of sparse-effect loci (nγ) for IC1–IC3 and selected variables. Estimates for each of the three MCMC chains are in Table S10.",
  "Table S10" = "BSLMM estimates of PVE, PGE and number of sparse-effect loci (nγ) for each of the three MCMC chains, for the variables and sets in Table S9.",
  "Table S11" = "Associated markers and their priority score for candidate-gene search: 2 if detected within a species, plus 1 each for a top WZA window, association with more than one variable and detection by LFMM (possible range 0–5; observed maximum 3), with the evidence behind each score.",
  "Table S12" = "Gene and mRNA features in the TDr96 annotation within 10 kb of markers scoring at least 2 in Table S11, with the associated variables, analysis set and distance from the marker (0, marker within the gene). For the two multi-variable full-panel markers annotated in detail, the position class within the gene and the best Araport11 (Arabidopsis thaliana) match from reciprocal BLASTP are also given.",
  "Table S13" = "Per-accession genomic offsets (mean and among-model standard deviation across six CMIP6 models, SSP5-8.5, 2081–2100) for the full panel and within D. rotundata and D. praehensilis, with Köppen–Geiger classes for 1991–2020 and 2071–2099 at each collection site.",
  "Table S14" = "Accessions in the two 50-accession cores compared as training sets: the core selected from genotypes with Core Hunter and used as the training set (23 D. rotundata, 8 western and 7 Cameroonian D. praehensilis, 3 D. abyssinica and 9 hybrids) and the original core (18 D. rotundata, 13 western D. praehensilis, 12 D. abyssinica and 7 hybrids; no Cameroonian D. praehensilis), with genetic group and collection coordinates. Eleven accessions are in both cores.",
  "Table S15" = "Prediction accuracy of the original and genotype-only 50-accession cores used as training sets.",
  "Table S16" = "Within-group genomic estimated adaptive values (GEAVs; out-of-fold predictions expressed as deviations from the genetic-group mean) for 167 accessions, 36 environmental variables and three ICs (Table S17 lists the top-ranked candidate parents).",
  "Table S17" = "Top-ranked candidate parents within each genetic group for low soil pH (PHIHOX) and high annual precipitation (BIO12), with their within-group GEAVs (Table S16).",
  "Table S18" = "Environmental Genome Wide Association Significant Markers with allelic effect (β), standard error, allele frequency, proportion of variance explained (PVE), q-value and genomic inflation factor (λ), for the full panel and within D. rotundata and D. praehensilis.",
  "Table S19" = "GenoBridge best model and hold-out accuracy (34 accessions) for all 39 variables, compared with the best model and the structure-only baseline under random 10-fold cross-validation.",
  NULL)
captions <- captions[!vapply(captions, is.null, logical(1))]
build <- list(
  "Table S1"  = function() { v <- rdt("TableS1_environmental_variables_VIF.csv"); d <- rdt("env_variable_definitions.csv")
                  if (is.null(v)) return(d); names(v)[1] <- "Variable"
                  v <- data.frame(Variable = v$Variable, `Retained after VIF pruning (VIF < 10)` = v[[2]],
                                  `Used in genomic offset` = if (ncol(v) > 2) v[[3]] else NA, check.names = FALSE)
                  if (is.null(d)) v else merge(d, v, by = "Variable", all = TRUE, sort = FALSE) },
  "Table S2"  = function() rdt("TableS2_missing_filter_by_group.csv"),
  "Table S3"  = function() rdt(file.path("tables", "comparison", "TableS3_strict_vs_relaxed.csv")),
  "Table S4"  = function() rdt("TableS4_IBD_IBE_MMRR_partialMantel.csv"),
  "Table S5"  = function() rdt("TableS5_IC_loadings.csv"),
  "Table S6"  = function() rdt("TableS6_WZA_top_windows.csv"),
  "Table S7"  = function() rdt("TableS7_lambda_and_counts.csv"),
  "Table S8"  = function() rdt("TableS8_LFMM_significant.csv"),
  "Table S9"  = function() rdt("TableS9_BSLMM.csv"),
  "Table S10" = function() rdt("TableS10_BSLMM_by_chain.csv"),
  "Table S11" = function() rdt("TableS11_candidate_markers_prioritised.csv"),
  "Table S12" = function() rdt("TableS12_candidate_genes.csv"),
  "Table S13" = function() { a <- rdt("TableS13_offset_per_accession.csv"); k <- rdt("TableS13b_offset_with_Koppen.csv")
                  if (!is.null(a) && !is.null(k)) a <- merge(a, k[, intersect(c("SampleID", "KG_1991_2020", "KG_2071_2099", "changed"), names(k))],
                                                            by = "SampleID", all.x = TRUE, sort = FALSE)
                  with_names(a) },
  "Table S14" = function() { fs <- list.files(TAB, "^TableS14_core_.*\\.csv$", full.names = TRUE)
                  fs <- fs[!grepl("collection_selected", fs)]; if (!length(fs)) return(with_names(rdt("TableS14_core_collection_selected.csv")))
                  sel <- tryCatch(rdt("TableS15_core_comparison.csv"), error = function(e) NULL)
                  best <- if (!is.null(sel) && "selected" %in% names(sel)) sel$core[sel$selected %in% c(TRUE, "TRUE")][1] else NA
                  cores <- lapply(fs, read.csv); names(cores) <- sub("^TableS14_core_(.*)\\.csv$", "\\1", basename(fs))
                  ids <- unique(unlist(lapply(cores, `[[`, "SampleID")))
                  out <- data.frame(SampleID = ids, `Genetic group` = meta_all$Population.assignment[match(ids, meta_all$SampleID)],
                                    Longitude = meta_all$Longitude[match(ids, meta_all$SampleID)],
                                    Latitude = meta_all$Latitude[match(ids, meta_all$SampleID)], check.names = FALSE)
                  for (k in names(cores)) out[[paste0(ifelse(identical(k, best), "Selected core: ", "Core: "), gsub("_", "-", k))]] <-
                    ifelse(ids %in% cores[[k]]$SampleID, "yes", "")
                  with_names(out) },
  "Table S15" = function() { d <- rdt("TableS15_core_comparison_detail.csv"); if (is.null(d)) rdt("TableS15_core_comparison.csv") else d },
  "Table S16" = function() with_names(rdt("TableS16_GEAV_within_group.csv")),
  "Table S17" = function() with_names(rdt("TableS17_candidate_parents_within_group.csv")),
  "Table S18" = function() rdt("TableS18_EGWAS_significant_with_effects.csv"),
  "Table S19" = function() rdt("TableS19_GenoBridge.csv"))
wb <- openxlsx::createWorkbook(); hs <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9E2F3")
index <- data.frame(Sheet = character(), Caption = character(), Rows = integer())
for (sh in names(build)) {
  x <- tryCatch(build[[sh]](), error = function(e) { msg("Problem building", sh, "-", conditionMessage(e)); NULL })
  if (is.null(x)) { msg("WARNING:", sh, "has no data - check the script that writes it"); x <- data.frame(Note = "not yet produced") }
  openxlsx::addWorksheet(wb, sh); openxlsx::writeData(wb, sh, x, headerStyle = hs); openxlsx::freezePane(wb, sh, firstRow = TRUE)
  index[nrow(index) + 1, ] <- list(sh, if (!is.null(captions[[sh]])) captions[[sh]] else "", nrow(x))
}
openxlsx::addWorksheet(wb, "Index"); openxlsx::writeData(wb, "Index", index, headerStyle = hs)
openxlsx::setColWidths(wb, "Index", 2, 120); openxlsx::worksheetOrder(wb) <- c(length(build) + 1, seq_along(build))
openxlsx::saveWorkbook(wb, tab("Supplementary_Tables_S1-S19.xlsx"), overwrite = TRUE)
msg("Wrote", tab("Supplementary_Tables_S1-S19.xlsx")); print(index[, c("Sheet", "Rows")])

rd <- function(f) { p <- tab(f); x <- if (file.exists(p)) safe_read(p) else NULL; if (!is.null(x) && nrow(x)) x else NULL }
L <- c("# Values for manuscript placeholders", paste("Genotype set:", CFG$geno_set), "")
add <- function(...) L <<- c(L, paste0(...))
if (!is.null(x <- rd("Table1_groups_He_ranges.csv"))) { add("## Results - Panel (Table 1)")
  add("- He range: ", paste(round(range(x$He), 3), collapse = " to ")); add("- per group:\n", paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("TableS2_missing_filter_by_group.csv"))) add("- Missing-filter loss by group (%SNPs with missing): ", paste(x$group, round(x$pct_SNPs_with_missing, 1), collapse = "; "))
if (!is.null(x <- rd("FigS4_envPCA_variance.csv"))) add("- Env PCA: PC1 ", x$var_pct[1], "%, PC2 ", x$var_pct[2], "%")
if (!is.null(x <- rd("TableS4_IBD_IBE_MMRR_partialMantel.csv"))) { add("## IBD / IBE (MMRR)"); add(paste(capture.output(print(x[, c("set","term","coef","p","R2")])), collapse = "\n")) }
if (!is.null(x <- rd("Table2_pRDA_summary.csv"))) { add("## Partial RDA (Table 2)"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("FST_pairwise_genomewide.csv"))) add("## FST pairwise: ", paste(x$contrast, round(x$WC_FST, 3), collapse = "; "))
if (!is.null(x <- rd("EGWAS_summary_by_set.csv"))) { add("## EGWAS"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("TableS6_WZA_top_windows.csv"))) add("- WZA top windows: ", paste(names(table(x$set)), table(x$set), collapse = "; "), "; overlapping EGWAS: ", sum(x$overlaps_EGWAS))
if (!is.null(x <- rd("TableS8_LFMM_significant.csv"))) { k <- tapply(x$SNP, x$set, function(v) length(unique(v)))
  add("- LFMM significant SNPs: ", paste(names(k), k, sep = " = ", collapse = "; ")) }
if (!is.null(x <- rd("TableS9_BSLMM.csv"))) { add("## BSLMM"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("FigS7_power.csv"))) { add("## Power (PVE at which power >= 0.8)")
  add(paste(sapply(split(x, x$set), function(d) paste0(d$set[1], ": ", ifelse(any(d$power >= .8), 100 * min(d$PVE[d$power >= .8]), ">20"), "%")), collapse = "; ")) }
if (!is.null(x <- rd("TableS13_offset_per_accession.csv"))) add("## Offset per accession: range ", paste(signif(range(x$offset_mean), 3), collapse = " to "))
if (!is.null(x <- rd("Koppen_transitions_sites.csv"))) add("- Sites changing Köppen class: ", sum(x$n_sites[x$changed]), " of ", sum(x$n_sites))
if (!is.null(x <- rd("FigS9_prediction_accuracy.csv"))) { add("## Prediction (mean r and gain over baseline)")
  add(paste(capture.output(print(x %>% group_by(scheme, model) %>% summarise(r = round(mean(r_pooled), 2), gain = round(mean(gain_over_baseline), 2), .groups = "drop"))), collapse = "\n")) }
if (!is.null(x <- rd("TableS17_candidate_parents_within_group.csv"))) { add("## Candidate parents"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("Fig2A_private_alleles_by_group.csv"))) { add("## Private alleles (raw / rarefied)"); add(paste(x$group, x$private_raw, round(x$private_rarefied, 1), collapse = "; ")) }
if (!is.null(x <- rd("Fig2B_froh_by_individual.csv"))) add("- FROH (ROH, corrected): max ", signif(max(x$FROH), 3), "; group means ",
    paste(tapply(signif(x$FROH, 3), x$group, mean), collapse = " / "))
if (!is.null(x <- rd("LD_decay_summary.csv"))) { add("## LD decay"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("FigS17_marker_ids.csv"))) add("- Fig. S17 markers: ", paste(x$SNP, collapse = ", "))
if (!is.null(x <- rd("TableS12_candidate_genes.csv"))) add("- Candidate genes found: ", nrow(x), " features near ", length(unique(x$SNP)), " markers")
writeLines(L, tab("placeholder_values.md")); msg("Wrote", tab("placeholder_values.md"))
