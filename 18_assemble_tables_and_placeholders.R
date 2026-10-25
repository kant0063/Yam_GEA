################################################################################
# 18_assemble_tables_and_placeholders.R                    Analysis Plan: A10
# Collects every table into Supplementary_Tables_revision.xlsx and writes
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
files <- list.files(TAB, pattern = "\\.csv$", full.names = TRUE)
# consistency check: Fig. S14b / Table S3 must come from the 13b run that also wrote Table S3b
f3b <- tab("TableS3b_core_comparison.csv"); f14b <- tab("FigS14b_core_vs_random_training.csv")
if (file.exists(f3b) && file.exists(f14b) && file.mtime(f14b) + 60 < file.mtime(f3b))
  warning("FigS14b_core_vs_random_training.csv is older than TableS3b - rerun 13b_core_comparison.R before 18")
wb <- openxlsx::createWorkbook(); used <- character(0)
index <- data.frame(sheet = character(), file = character(), rows = integer(), cols = integer(), note = character())
for (f in files) {
  x <- safe_read(f)
  base <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", sub("\\.csv$", "", basename(f)))
  sh <- substr(base, 1, 31); k <- 1
  while (tolower(sh) %in% tolower(used)) { k <- k + 1; sh <- paste0(substr(base, 1, 31 - nchar(k) - 1), "_", k) }  # unique, <= 31 chars
  if (is.null(x)) { index[nrow(index) + 1, ] <- list(NA, basename(f), 0L, 0L, "empty - no rows (e.g. no significant results)"); next }
  openxlsx::addWorksheet(wb, sh); openxlsx::writeData(wb, sh, x); used <- c(used, sh)
  index[nrow(index) + 1, ] <- list(sh, basename(f), nrow(x), ncol(x), "")
}
openxlsx::addWorksheet(wb, "INDEX"); openxlsx::writeData(wb, "INDEX", index)
openxlsx::worksheetOrder(wb) <- c(length(used) + 1, seq_along(used))
print(index)
openxlsx::saveWorkbook(wb, tab("Supplementary_Tables_revision.xlsx"), overwrite = TRUE)

rd <- function(f) { p <- tab(f); x <- if (file.exists(p)) safe_read(p) else NULL; if (!is.null(x) && nrow(x)) x else NULL }
L <- c("# Values for manuscript placeholders", paste("Genotype set:", CFG$geno_set), "")
add <- function(...) L <<- c(L, paste0(...))
if (!is.null(x <- rd("Table1_groups_He_ranges.csv"))) { add("## Results - Panel (Table 1)")
  add("- He range: ", paste(round(range(x$He), 3), collapse = " to ")); add("- per group:\n", paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("TableS7a_missing_filter_by_group.csv"))) add("- Missing-filter loss by group (%SNPs with missing): ", paste(x$group, round(x$pct_SNPs_with_missing, 1), collapse = "; "))
if (!is.null(x <- rd("FigS17_envPCA_variance.csv"))) add("- Env PCA: PC1 ", x$var_pct[1], "%, PC2 ", x$var_pct[2], "%")
if (!is.null(x <- rd("TableS8_IBD_IBE_MMRR_partialMantel.csv"))) { add("## IBD / IBE (MMRR)"); add(paste(capture.output(print(x[, c("set","term","coef","p","R2")])), collapse = "\n")) }
if (!is.null(x <- rd("Table2_pRDA_summary.csv"))) { add("## Partial RDA (Table 2)"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("FST_pairwise_genomewide.csv"))) add("## FST pairwise: ", paste(x$contrast, round(x$WC_FST, 3), collapse = "; "))
if (!is.null(x <- rd("EGWAS_summary_by_set.csv"))) { add("## EGWAS"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("TableS9_WZA_top_windows.csv"))) add("- WZA top windows: ", paste(names(table(x$set)), table(x$set), collapse = "; "), "; overlapping EGWAS: ", sum(x$overlaps_EGWAS))
if (!is.null(x <- rd("TableS10_LFMM_significant.csv"))) { k <- tapply(x$SNP, x$set, function(v) length(unique(v)))
  add("- LFMM significant SNPs: ", paste(names(k), k, sep = " = ", collapse = "; ")) }
if (!is.null(x <- rd("TableS11_BSLMM.csv"))) { add("## BSLMM"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("FigS18_power.csv"))) { add("## Power (PVE at which power >= 0.8)")
  add(paste(sapply(split(x, x$set), function(d) paste0(d$set[1], ": ", ifelse(any(d$power >= .8), 100 * min(d$PVE[d$power >= .8]), ">20"), "%")), collapse = "; ")) }
if (!is.null(x <- rd("TableS6_offset_per_accession.csv"))) add("## Offset per accession: range ", paste(signif(range(x$offset_mean), 3), collapse = " to "))
if (!is.null(x <- rd("Koppen_transitions_sites.csv"))) add("- Sites changing Köppen class: ", sum(x$n_sites[x$changed]), " of ", sum(x$n_sites))
if (!is.null(x <- rd("FigS14_prediction_accuracy.csv"))) { add("## Prediction (mean r and gain over baseline)")
  add(paste(capture.output(print(x %>% group_by(scheme, model) %>% summarise(r = round(mean(r_pooled), 2), gain = round(mean(gain_over_baseline), 2), .groups = "drop"))), collapse = "\n")) }
if (!is.null(x <- rd("candidate_parents_within_group.csv"))) { add("## Candidate parents"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("Fig2A_private_alleles_by_group.csv"))) { add("## Private alleles (raw / rarefied)"); add(paste(x$group, x$private_raw, round(x$private_rarefied, 1), collapse = "; ")) }
if (!is.null(x <- rd("Fig2B_froh_by_individual.csv"))) add("- FROH (ROH, corrected): max ", signif(max(x$FROH), 3), "; group means ",
    paste(tapply(signif(x$FROH, 3), x$group, mean), collapse = " / "))
if (!is.null(x <- rd("LD_decay_summary.csv"))) { add("## LD decay"); add(paste(capture.output(print(x)), collapse = "\n")) }
if (!is.null(x <- rd("FigS13_marker_ids.csv"))) add("- Fig. S13 markers: ", paste(x$SNP, collapse = ", "))
if (!is.null(x <- rd("candidate_genes.csv"))) add("- Candidate genes found: ", nrow(x), " features near ", length(unique(x$SNP)), " markers")
writeLines(L, tab("placeholder_values.md")); msg("Wrote", tab("placeholder_values.md"))
