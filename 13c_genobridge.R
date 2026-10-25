################################################################################
# 13c_genobridge.R                                         Analysis Plan: A9 (GenoBridge)
# GenoBridge v1.0.1 (Singh & Mendu 2026, Plant Methods; github.com/nps-genomics/GenoBridge)
# runs from the command line on LINUX x86_64 only (conda package). It selects among
# ridge, random forest, XGBoost, neural-network and stacking models for each trait and
# uses its OWN internal hold-out validation, so it cannot be run on our blocked folds;
# its accuracy is compared with our random cross-validation.
#
# Step 1 (here, any OS):  Rscript 13c_genobridge.R      -> writes the input files + commands
# Step 2 (Linux / HPC / Docker): run the printed genobridge-predict command
# Step 3 (here):          Rscript 13c_genobridge.R      -> detects the results and writes
#                         TableS14c_GenoBridge.csv (compared with Fig. S14 random CV)
################################################################################
source("00_config_and_helpers.R"); M <- load_master()
gbdir <- res("genobridge"); dir.create(gbdir, showWarnings = FALSE)
out_dir <- file.path(gbdir, "results")
IC <- readRDS(res("ICs.rds")); Y <- scale(cbind(M$E, IC)); ids <- rownames(Y)
ph <- data.frame(accession_id = ids, Y, check.names = FALSE)
write.csv(ph, file.path(gbdir, "phenotypes.csv"), row.names = FALSE, na = "")
vcf <- normalizePath(CFG$vcf_strict, mustWork = FALSE)
msg("Phenotype file:", file.path(gbdir, "phenotypes.csv"), "|", ncol(Y), "variables,", nrow(Y), "accessions")
msg("VCF sample IDs must equal accession_id (the same IDs as in the master object).")
cat("\n# --- Linux (conda) -------------------------------------------------------------\n",
    "conda create -n genobridge --override-channels -c nps-genomics -c conda-forge -c bioconda python=3.12 genobridge=1.0.1 --solver=libmamba -y\n",
    "conda activate genobridge\n",
    sprintf("genobridge-predict --phenotype %s --vcf %s --accession-col accession_id --no-outlier-removal --output %s\n",
            shQuote(normalizePath(file.path(gbdir, "phenotypes.csv"))), shQuote(vcf), shQuote(normalizePath(out_dir, mustWork = FALSE))),
    "\n# --- macOS: run the same commands inside a Linux container -------------------------\n",
    sprintf("docker run --platform linux/amd64 -it -v %s:/work -w /work continuumio/miniconda3 bash\n",
            shQuote(normalizePath(getwd()))),
    "# (inside the container, use paths relative to /work)\n", sep = "")

# ---- parse results once they exist ------------------------------------------------
rf <- file.path(out_dir, "phenotype_prediction_results.csv")
if (file.exists(rf)) {
  gb <- read.csv(rf, check.names = FALSE, stringsAsFactors = FALSE)
  msg("GenoBridge results columns:", paste(names(gb), collapse = ", "))
  col_tr <- grep("trait|phenotype", names(gb), ignore.case = TRUE, value = TRUE)[1]
  col_r  <- grep("^best_?r$|best.*corr|best.*r$", names(gb), ignore.case = TRUE, value = TRUE)[1]
  col_m  <- grep("best.*model", names(gb), ignore.case = TRUE, value = TRUE)[1]
  if (any(is.na(c(col_tr, col_r)))) stop("Could not find trait / best-r columns - check the column names printed above")
  gbt <- data.frame(trait = gb[[col_tr]], GenoBridge_r = as.numeric(gb[[col_r]]),
                    GenoBridge_model = if (!is.na(col_m)) gb[[col_m]] else NA)
  acc <- read.csv(tab("FigS14_prediction_accuracy.csv"))
  ours <- acc %>% filter(scheme == "random 10-fold") %>% group_by(trait) %>%
    summarise(best_ours_r = max(r_pooled[model != "baseline"], na.rm = TRUE),
              best_ours_model = model[model != "baseline"][which.max(r_pooled[model != "baseline"])],
              baseline_r = r_pooled[model == "baseline"][1], .groups = "drop")
  tabc <- left_join(gbt, ours, by = "trait") %>% mutate(diff_vs_ours = GenoBridge_r - best_ours_r,
                                                       gain_over_baseline = GenoBridge_r - baseline_r)
  write.csv(tabc, tab("TableS14c_GenoBridge.csv"), row.names = FALSE)
  print(tabc %>% summarise(n = n(), GenoBridge_mean_r = mean(GenoBridge_r, na.rm = TRUE),
                           ours_mean_r = mean(best_ours_r, na.rm = TRUE), mean_diff = mean(diff_vs_ours, na.rm = TRUE),
                           mean_gain_over_baseline = mean(gain_over_baseline, na.rm = TRUE)))
  if (!all(is.na(gbt$GenoBridge_model))) print(table(gbt$GenoBridge_model))
} else msg("No GenoBridge results yet in", out_dir, "- run the command above, then rerun this script.")
