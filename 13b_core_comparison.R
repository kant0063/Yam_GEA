################################################################################
# 13b_core_comparison.R                                    Analysis Plan: A9 (core)
# Compares candidate 50-accession cores as training sets and keeps the one with the
# highest prediction accuracy (reported in the manuscript):
#   - "original": the Table S14 core (CFG$core_csv; a title line above the header is handled)
#   - "genotype_only": a core re-selected from genotypes only with Core Hunter
#     (needs the corehunter R package + Java; skipped if unavailable)
# Each core is benchmarked as in script 13: train on the core, predict the non-core
# accessions, compared with random and group-stratified training sets of the same size.
# Reuses the model set-up of script 13 (no need to rerun 13).
# Outputs: TableS15_core_comparison.csv, TableS14_core_collection_selected.csv (best core),
#          FigS9D_core_vs_random_training.csv/.png (best core).
################################################################################
src <- readLines("13_genomic_prediction_revised.R")
a <- grep('^source\\("00_config_and_helpers.R"\\)', src)[1]; b <- grep("^run_cv <- function", src)[1] - 1
eval(parse(text = src[a:b]))                       # data, kernels, folds and fit_pred() from script 13
grp <- meta$Population.assignment[match(ids, meta$SampleID)]
models_light <- c("baseline", "rrBLUP", "GBLUP_gauss", "GBLUP_exp")
if (is.null(CFG$core_reps)) CFG$core_reps <- 50
if (is.null(CFG$core_csv))  CFG$core_csv  <- "Yam_core_collection.csv"

read_core <- function(f) {                          # tolerate a title line above the header
  L <- readLines(f, warn = FALSE); L <- L[nzchar(trimws(gsub(",", "", L)))]
  v <- trimws(gsub('"', "", sapply(strsplit(L, ","), `[`, 1)))
  intersect(v, ids)
}
cores <- list()
if (file.exists(CFG$core_csv)) cores$original <- read_core(CFG$core_csv)
if (requireNamespace("corehunter", quietly = TRUE)) {
  set.seed(CFG$seed)
  gd <- corehunter::genotypes(G, format = "biparental")
  ch <- corehunter::sampleCore(corehunter::coreHunterData(genotypes = gd), size = 50,
          obj = list(corehunter::objective("EN", "MR"), corehunter::objective("AN", "MR")))
  cores$genotype_only <- intersect(ch$sel, ids)
} else msg("corehunter not installed - only the original core is evaluated")
stopifnot(length(cores) > 0)
for (k in names(cores)) msg("Core", k, ":", length(cores[[k]]), "accessions; groups:",
                            paste(names(table(grp[match(cores[[k]], ids)])), table(grp[match(cores[[k]], ids)]), collapse = ", "))

acc_sub <- function(pred, y, g) c(r_pooled = cor(pred, y),
  r_within_group = mean(sapply(split(seq_along(y), g), function(i) if (length(i) > 4) cor(pred[i], y[i]) else NA), na.rm = TRUE))
eval_core <- function(core_ids) {
  core_i <- match(core_ids, ids); nc <- length(core_i); comp <- table(grp[core_i])
  strat <- function() unlist(lapply(names(comp), function(g) { pool <- which(grp == g)
    pool[sample.int(length(pool), min(comp[[g]], length(pool)))] }))
  set.seed(CFG$seed)
  sets <- list(core = list(core_i),
               random_matched = replicate(CFG$core_reps, sample.int(length(ids), nc), simplify = FALSE),
               stratified_matched = replicate(CFG$core_reps, strat(), simplify = FALSE))
  out <- list()
  for (trait in colnames(Y)) for (mo in models_light) for (sch in names(sets)) {
    st <- sapply(sets[[sch]], function(tr) { te <- setdiff(seq_along(ids), tr)
      acc_sub(fit_pred(mo, Y[, trait], tr, te), Y[te, trait], grp[te]) })
    out[[paste(trait, mo, sch)]] <- data.frame(trait, model = mo, training = sch, n_train = nc,
      r_pooled = mean(st["r_pooled", ]), r_pooled_sd = if (ncol(st) > 1) sd(st["r_pooled", ]) else NA,
      r_within_group = mean(st["r_within_group", ], na.rm = TRUE))
  }
  ca <- bind_rows(out)
  base <- ca %>% filter(model == "baseline") %>% dplyr::select(trait, training, r_base = r_pooled)
  ca %>% left_join(base, by = c("trait", "training")) %>% mutate(gain_over_baseline = r_pooled - r_base)
}
res_all <- lapply(names(cores), function(k) { msg("Evaluating core:", k); eval_core(cores[[k]]) %>% mutate(core = k) })
names(res_all) <- names(cores)

cmp <- bind_rows(res_all) %>% filter(training == "core", model != "baseline") %>%
  group_by(core) %>% summarise(n = first(n_train), mean_r = mean(r_pooled, na.rm = TRUE),
                               mean_r_within = mean(r_within_group, na.rm = TRUE),
                               mean_gain = mean(gain_over_baseline, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_r))
best <- cmp$core[1]
cmp$selected <- cmp$core == best
# Table S15: summary of every core (composition + accuracy by training design and model)
comp <- sapply(names(cores), function(k) paste(names(table(grp[match(cores[[k]], ids)])), table(grp[match(cores[[k]], ids)]), collapse = "; "))
s3b <- bind_rows(res_all) %>% group_by(core, training, model) %>%
  summarise(n_train = first(n_train), r = mean(r_pooled, na.rm = TRUE), r_within = mean(r_within_group, na.rm = TRUE),
            gain_over_baseline = mean(gain_over_baseline, na.rm = TRUE), .groups = "drop") %>%
  mutate(composition = comp[core], selected = core == best)
write.csv(cmp, tab("TableS15_core_comparison.csv"), row.names = FALSE)
write.csv(s3b, tab("TableS15_core_comparison_detail.csv"), row.names = FALSE)
for (k in names(cores)) write.csv(meta[match(cores[[k]], meta$SampleID), c("SampleID", "Population.assignment", "Longitude", "Latitude")],
                                  tab(paste0("TableS14_core_", k, ".csv")), row.names = FALSE)
msg("Core with the highest prediction accuracy:", best); print(cmp)

ca <- res_all[[best]]
write.csv(ca, tab("FigS9D_core_vs_random_training.csv"), row.names = FALSE)
write.csv(meta[match(cores[[best]], meta$SampleID), c("SampleID", "Population.assignment", "Longitude", "Latitude")],
          tab("TableS14_core_collection_selected.csv"), row.names = FALSE)
print(ca %>% group_by(training, model) %>% summarise(r = round(mean(r_pooled), 2), r_within = round(mean(r_within_group, na.rm = TRUE), 2),
                                                    gain = round(mean(gain_over_baseline), 2), .groups = "drop"), n = Inf)
z <- ca %>% filter(training %in% c("core", "random_matched")) %>% dplyr::select(trait, model, training, r_pooled, r_pooled_sd) %>%
  tidyr::pivot_wider(names_from = training, values_from = c(r_pooled, r_pooled_sd))
print(z %>% group_by(model) %>% summarise(core_beats_random = mean(r_pooled_core > r_pooled_random_matched, na.rm = TRUE),
  mean_z = round(mean((r_pooled_core - r_pooled_random_matched) / r_pooled_sd_random_matched, na.rm = TRUE), 1), .groups = "drop"))
lab_tr <- c(core = "Core", random_matched = "Random\n(size-matched)", stratified_matched = "Group-stratified\n(size-matched)")
ca$training <- factor(ca$training, levels = names(lab_tr), labels = lab_tr)
p <- ggplot(filter(ca, model != "baseline"), aes(training, r_pooled, fill = training)) + geom_boxplot(outlier.size = .6) +
  facet_wrap(~model) + geom_hline(yintercept = 0, linetype = 2) +
  labs(x = NULL, y = "Accuracy on non-training accessions (r)",
       title = NULL) +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(size = 8))
ggsave(fig("FigS9D_core_vs_random_training.png"), p, width = 9, height = 4.5, dpi = 300, bg = "white")
