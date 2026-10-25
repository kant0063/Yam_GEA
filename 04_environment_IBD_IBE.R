################################################################################
# 04_environment_IBD_IBE.R                                 Analysis Plan: A3 + A4
# Fig. S4 (environmental PCA coloured by genetic group + site map over PC1),
# Table S4 (isolation by distance vs environment: MMRR + partial Mantel), for the
# full panel and within D. rotundata / D. praehensilis.
#
# v2 fixes: "subscript out of bounds" came from NA sample IDs returned by
# subsets() for accessions whose Population.assignment names no species
# (Species_group = NA). subsets() in 00 now uses %in%; this script also checks
# IDs, variables and the thinning index before each analysis, drops
# zero-variance variables within a subset, labels unassigned samples in the
# plots (the 14 "removed rows"), and uses world_borders() (ggplot2 >= 4.0).
################################################################################
source("00_config_and_helpers.R"); need(c("vegan", "geosphere", "maps")); M <- load_master(); meta <- M$meta

# ---- checks ------------------------------------------------------------------
meta$Species_group[is.na(meta$Species_group)] <- "unassigned"
if (any(meta$Species_group == "unassigned")) {
  msg("Samples with no species in Population.assignment:", sum(meta$Species_group == "unassigned"))
  print(table(meta$Population.assignment[meta$Species_group == "unassigned"], useNA = "ifany"))
}
vars <- intersect(M$vif_keep, colnames(M$E))
if (length(vars) < length(M$vif_keep)) msg("WARNING: VIF variables missing from M$E:", paste(setdiff(M$vif_keep, colnames(M$E)), collapse = ", "))
if (length(vars) < 2) stop("Fewer than 2 VIF-retained variables found in M$E - rerun 01")
thin <- M$thin
if (length(thin) != ncol(M$G)) { msg("Thinning index length does not match genotypes - recomputing")
  thin <- thin_snps(M$snp[match(colnames(M$G), M$snp$SNP), ]) }
bad_xy <- is.na(meta$Longitude) | is.na(meta$Latitude)
if (any(bad_xy)) msg("WARNING:", sum(bad_xy), "samples without coordinates are excluded from IBD/IBE")

# ---- A3: environmental PCA and map ------------------------------------------
Es <- scale(M$E[, apply(M$E, 2, sd, na.rm = TRUE) > 0]); Es[is.na(Es)] <- 0
epca <- prcomp(Es); ev <- round(100 * epca$sdev^2 / sum(epca$sdev^2), 1)
d3 <- data.frame(meta, PC1 = epca$x[match(meta$SampleID, rownames(Es)), 1], PC2 = epca$x[match(meta$SampleID, rownames(Es)), 2])
load <- data.frame(var = rownames(epca$rotation), epca$rotation[, 1:2] * 4)
load_top <- load[order(-(load$PC1^2 + load$PC2^2)), ][1:min(15, nrow(load)), ]   # label the 15 longest loadings only
pA <- ggplot(d3, aes(PC1, PC2, colour = Population.assignment)) + geom_point(size = 2, alpha = .8) +
  geom_segment(data = load, aes(0, 0, xend = PC1, yend = PC2), inherit.aes = FALSE, colour = "grey50", linewidth = .3) +
  { if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(data = load_top, aes(PC1, PC2, label = var), inherit.aes = FALSE, size = 2.6, colour = "grey20",
                               max.overlaps = Inf, segment.colour = "grey60", min.segment.length = 0)
    else geom_text(data = load_top, aes(PC1, PC2, label = var), inherit.aes = FALSE, size = 2.6, colour = "grey20", check_overlap = TRUE) } +
  labs(x = paste0("Env PC1 (", ev[1], "%)"), y = paste0("Env PC2 (", ev[2], "%)"), colour = "Group") + theme_bw()
xl <- range(d3$Longitude, na.rm = TRUE); yl <- range(d3$Latitude, na.rm = TRUE)
shapes <- c(rotundata = 21, praehensilis = 22, abyssinica = 24, hybrid = 23, unassigned = 25)
pB <- ggplot(d3[!bad_xy, ], aes(Longitude, Latitude)) +
  world_borders(xlim = xl + c(-2, 2), ylim = yl + c(-2, 2), colour = "grey60") +
  geom_point(aes(fill = PC1, shape = Species_group), size = 2.6) +
  scale_shape_manual(values = shapes[names(shapes) %in% unique(d3$Species_group)]) +
  scale_fill_viridis_c(name = "Env PC1") + coord_quickmap(xlim = xl + c(-1, 1), ylim = yl + c(-1, 1)) +
  labs(x = "Longitude", y = "Latitude", shape = "Group") + theme_bw()
save_panel(pA, "FigS4A"); save_panel(pB, "FigS4B")   # composite Figure S4 is built by script 17
ggsave(fig("FigS4A_environmental_PCA.png"), pA, width = 8, height = 6, dpi = 300)
ggsave(fig("FigS4B_site_map_envPC1.png"), pB, width = 7, height = 5, dpi = 300)
write.csv(data.frame(PC = 1:5, var_pct = ev[1:5]), tab("FigS4_envPCA_variance.csv"), row.names = FALSE)
write.csv(d3 %>% group_by(Population.assignment) %>% summarise(n = n(), envPC1_mean = mean(PC1), envPC1_min = min(PC1),
          envPC1_max = max(PC1)), tab("FigS4_envPC1_by_group.csv"), row.names = FALSE)

# ---- A4: IBD vs IBE ---------------------------------------------------------
run_ibe <- function(ids, label) {
  ids <- ids[!is.na(ids)]
  ids <- Reduce(intersect, list(ids, rownames(M$E), rownames(M$G), meta$SampleID[!bad_xy]))
  if (length(ids) < 10) { msg("Skipping", label, "- only", length(ids), "samples"); return(NULL) }
  m <- meta[match(ids, meta$SampleID), ]
  E <- M$E[ids, vars, drop = FALSE]; E <- E[, apply(E, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  msg(label, ":", length(ids), "samples,", ncol(E), "env variables,", sum(thin), "thinned SNPs")
  Dgeo <- geosphere::distm(as.matrix(m[, c("Longitude", "Latitude")])) / 1000
  Denv <- as.matrix(dist(scale(E)))
  Gs <- maf_filter(M$G[ids, thin, drop = FALSE]); Gs[is.na(Gs)] <- 0
  Dgen <- as.matrix(dist(scale(Gs, scale = FALSE))) / sqrt(ncol(Gs))
  r <- mmrr(Dgen, list(geography = Dgeo, environment = Denv))
  pm1 <- vegan::mantel.partial(as.dist(Dgen), as.dist(Denv), as.dist(Dgeo), permutations = CFG$perms)   # env | geo
  pm2 <- vegan::mantel.partial(as.dist(Dgen), as.dist(Dgeo), as.dist(Denv), permutations = CFG$perms)   # geo | env
  r$partial_mantel_r <- c(pm2$statistic, pm1$statistic); r$partial_mantel_p <- c(pm2$signif, pm1$signif)
  r$set <- label; r$n <- length(ids); r$n_sites <- n_distinct(m$SiteID); r
}
S <- subsets(M)
s8 <- bind_rows(lapply(names(S), function(s) run_ibe(S[[s]], s)))
write.csv(s8, tab("TableS4_IBD_IBE_MMRR_partialMantel.csv"), row.names = FALSE); print(s8)
