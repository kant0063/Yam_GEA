################################################################################
# 11_genomic_offset_revised.R                              Analysis Plan: A8
# REPLACES General_genetic_offset_script_with_regional_plots.R and panel B of
# Combined_RDA_GenomicOffset_Figure.R.
# Problems found in the originals (each fixed here):
#   1. Fig. 3B was the GRAND MEAN over SSP245 + SSP585 and 2041-2060 + 2081-2100,
#      although labelled SSP585 2081-2100 -> now one scenario/period (CFG$ssp, CFG$period).
#   2. Missing (ocean) climate values were mean-imputed, so ocean cells received
#      offsets -> NA cells are kept NA and masked.
#   3. Combined_RDA_GenomicOffset_Figure.R computed |future prediction of each grid
#      cell - MEAN current prediction across all samples|, i.e. distance from the
#      panel's average climate, not change at a location. This is why offsets were
#      lowest around the sampled sites. -> LEA::genetic.gap compares current vs
#      future climate AT THE SAME location.
#   4. No extrapolation control -> cells whose CURRENT climate is outside the sampled envelope
#      (ExDet NT1 < 0 or NT2 > 1) are masked; the share of cells with non-analog FUTURE climate
#      is reported (v2: masking on future novelty greyed out the entire map).
#   5. Colour scale: single sequential palette, low -> high stated in the legend.
# Outputs: Fig. 3B, per-GCM SD map, per-accession offsets (Table S13 columns),
#          within-species offsets.
################################################################################
source("00_config_and_helpers.R"); need(c("LEA", "terra", "geodata")); library(terra)
M <- load_master(); meta <- M$meta; S <- subsets(M)
vars <- M$vif_keep_bio                                 # bioclim variables only (future layers exist only for bio1-19)
bio_index <- as.integer(gsub("\\D", "", vars))
msg("Offset predictors:", paste(vars, collapse = ", "))

ext0 <- ext(min(meta$Longitude) - CFG$map_buffer, max(meta$Longitude) + CFG$map_buffer,
           min(meta$Latitude) - CFG$map_buffer, max(meta$Latitude) + CFG$map_buffer)
EV <- as.vector(ext0)
now <- crop(geodata::worldclim_global("bio", res = CFG$clim_res, path = CFG$clim_dir), ext0)[[bio_index]]
names(now) <- vars
pts <- vect(meta, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

# training environment = WorldClim at sample coordinates (same source as the future layers)
E_train <- as.matrix(terra::extract(now, pts)[, -1]); rownames(E_train) <- meta$SampleID
cells <- which(complete.cases(values(now)))           # land cells only - no imputation
E_now_grid <- values(now)[cells, , drop = FALSE]

f_lfmm <- res("geno_full.lfmm")
if (!file.exists(f_lfmm)) LEA::write.lfmm(maf_filter(M$G), f_lfmm)

grid_off <- acc_off <- fut_sites <- list(); novel <- rep(FALSE, length(cells)); fut_novel <- rep(FALSE, length(cells))
for (gcm in CFG$gcms) {
  msg("GCM", gcm)
  fut <- crop(geodata::cmip6_world(gcm, CFG$ssp, CFG$period, var = "bioc", res = CFG$clim_res, path = CFG$clim_dir), ext0)[[bio_index]]
  names(fut) <- vars
  E_fut_grid <- values(fut)[cells, , drop = FALSE]; ok <- complete.cases(E_fut_grid)
  E_fut_acc <- as.matrix(terra::extract(fut, pts)[, -1]); fut_sites[[gcm]] <- E_fut_acc
  go_grid <- rep(NA_real_, length(cells))
  go_grid[ok] <- LEA::genetic.gap(input = f_lfmm, env = E_train, new.env = E_now_grid[ok, ],
                                   pred.env = E_fut_grid[ok, ], K = CFG$K_lfmm_full, scale = TRUE)$offset
  grid_off[[gcm]] <- go_grid
  acc_off[[gcm]] <- LEA::genetic.gap(input = f_lfmm, env = E_train, new.env = E_train,
                                     pred.env = E_fut_acc, K = CFG$K_lfmm_full, scale = TRUE)$offset
  xd <- exdet(E_train, E_fut_grid); fut_novel <- fut_novel | !(xd$analog %in% TRUE)
}
novel <- !exdet(E_train, E_now_grid)$analog      # mask: CURRENT climate outside the sampled envelope
# Future climate under SSP5-8.5 2081-2100 is outside the sampled envelope almost everywhere (warming beyond the
# observed range), so masking on future novelty would grey out the whole map. Future novelty is reported instead.
msg("Land cells with non-analog FUTURE climate (any GCM):", round(100 * mean(fut_novel), 1), "%")
write.csv(data.frame(land_cells = length(cells), current_nonanalog_pct = 100 * mean(novel), future_nonanalog_pct = 100 * mean(fut_novel)),
          tab("offset_extrapolation_summary.csv"), row.names = FALSE)
# climate at sampling sites, current and future (GCM mean) -> Fig. S18 (script 17)
saveRDS(list(now = E_train, fut = Reduce(`+`, fut_sites) / length(fut_sites), group = meta$Population.assignment),
        res("climate_at_sites.rds"))
off_mean <- rowMeans(do.call(cbind, grid_off), na.rm = TRUE); off_sd <- apply(do.call(cbind, grid_off), 1, sd, na.rm = TRUE)

r_mean <- r_sd <- r_mask <- now[[1]]; values(r_mean) <- values(r_sd) <- values(r_mask) <- NA
r_mean[cells] <- ifelse(novel, NA, off_mean); r_sd[cells] <- ifelse(novel, NA, off_sd); r_mask[cells] <- as.numeric(novel)
writeRaster(c(r_mean, r_sd, r_mask), res(sprintf("offset_ssp%s_%s.tif", CFG$ssp, CFG$period)), overwrite = TRUE,
            names = c("offset_mean", "offset_sd_among_GCMs", "extrapolation_masked"))

# ---- Fig. 3B ---------------------------------------------------------------
df <- as.data.frame(c(r_mean, r_mask), xy = TRUE, na.rm = FALSE); names(df)[3:4] <- c("offset", "masked")
df <- df[!is.na(df$masked), ]
p <- ggplot() + geom_raster(data = df, aes(x, y, fill = offset)) +
  scale_fill_viridis_c(name = "Genomic offset\n(low \u2192 high)", na.value = "grey80") +
  world_borders(xlim = EV[1:2], ylim = EV[3:4], colour = "grey25", linewidth = .3) +
  geom_point(data = meta, aes(Longitude, Latitude), shape = 21, fill = "white", size = 1.6) +
  coord_quickmap(xlim = EV[1:2], ylim = EV[3:4], expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()
ggsave(fig(sprintf("Fig3B_offset_ssp%s_%s.png", CFG$ssp, CFG$period)), p, width = 8, height = 5.5, dpi = 300)

offset_summary <- data.frame(stat = c("grid_min", "grid_median", "grid_max", "grid_mean_sd_among_GCMs"),
  value = c(quantile(off_mean[!novel], c(0, .5, 1), na.rm = TRUE), mean(off_sd[!novel], na.rm = TRUE)))
write.csv(offset_summary, tab("offset_grid_summary.csv"), row.names = FALSE); print(offset_summary)
# ---- per-accession offsets (Table S13) and latitude trend ---------------------
acc <- data.frame(SampleID = meta$SampleID, group = meta$Population.assignment, Latitude = meta$Latitude,
                  offset_mean = rowMeans(do.call(cbind, acc_off)), offset_sd = apply(do.call(cbind, acc_off), 1, sd))
# within-species offsets (trained on each species only)
for (s in c("rotundata", "praehensilis")) {
  ids <- S[[s]]; fs <- res(paste0("geno_", s, ".lfmm")); if (!file.exists(fs)) LEA::write.lfmm(maf_filter(M$G[ids, ]), fs)
  o <- sapply(CFG$gcms, function(gcm) {
    fut <- crop(geodata::cmip6_world(gcm, CFG$ssp, CFG$period, var = "bioc", res = CFG$clim_res, path = CFG$clim_dir), ext0)[[bio_index]]
    Ef <- as.matrix(terra::extract(fut, pts)[, -1])[match(ids, meta$SampleID), ]
    LEA::genetic.gap(input = fs, env = E_train[ids, ], new.env = E_train[ids, ], pred.env = Ef,
                     K = CFG$K_lfmm_within[[s]], scale = TRUE)$offset })
  acc[[paste0("offset_within_", s)]] <- rowMeans(o)[match(acc$SampleID, ids)]
}
write.csv(acc, tab("TableS13_offset_per_accession.csv"), row.names = FALSE)
print(acc %>% group_by(group) %>% summarise(n = n(), median_offset = median(offset_mean), max_offset = max(offset_mean)))
print(summary(lm(offset_mean ~ Latitude, data = acc)))   # direction check for "increases northward"
msg("Grid cells masked as extrapolation:", sum(novel), "of", length(cells))
