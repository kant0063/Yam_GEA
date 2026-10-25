################################################################################
# 17_figures_main_and_supplementary.R                      Analysis Plan: A10
# Redraws every manuscript figure not already written by 02-16 from the revised
# results, so all figures come from the same data and settings:
#   Fig. 1A (PCA), Fig. 3, Fig. S15 (redrawn from Table S6), Fig. 4 (A partial RDA, B offset, C-D Köppen-Geiger),
#   Fig. S4, S5, S6, S7, S8, S9, S11, S13, S16.
# Figures produced elsewhere: 2A-B (14), 2C (06), S1 (02), S2-S3 (Yam_spatial_population_structure.R),
#   S10 (07), S12 (15), S14-S15 (13), S17 (04), S18 (10), 4A full/partial singles (05), 4B single (11).
# Each block is independent: if its inputs are missing it is skipped with a message.
# Needs: ggplot2; optional patchwork (Fig. 4 composite), igraph (Fig. S6), terra (Fig. 4B-D).
################################################################################
source("00_config_and_helpers.R"); M <- load_master(); meta <- M$meta
# On error, any figure file written in the failed block is removed (ggsave leaves a blank PNG if drawing fails)
block <- function(name, expr) {
  before <- list.files(FIG, full.names = TRUE); t0 <- Sys.time()
  tryCatch({ expr; msg("done:", name) }, error = function(e) {
    new <- setdiff(list.files(FIG, full.names = TRUE), before)
    new <- c(new, list.files(FIG, full.names = TRUE)[file.mtime(list.files(FIG, full.names = TRUE)) >= t0])
    if (length(new)) file.remove(unique(new))
    msg("SKIPPED", name, "-", conditionMessage(e)) })
}
vtype <- function(v) { k <- as.integer(gsub("\\D", "", v))
  ifelse(grepl("^IC", v), "climate IC", ifelse(grepl("^bio", v, ignore.case = TRUE),
         ifelse(!is.na(k) & k >= 12, "precipitation", "temperature"), "soil")) }
hits_all <- if (file.exists(tab("TableS4_EGWAS_significant_with_effects.csv"))) read.csv(tab("TableS4_EGWAS_significant_with_effects.csv")) else NULL
hits <- if (!is.null(hits_all)) filter(hits_all, set == "full") else NULL
Y <- if (file.exists(res("ICs.rds"))) cbind(M$E, readRDS(res("ICs.rds"))) else M$E

# ---------------- Fig. 1A ----------------
block("Fig1A", {
  d <- data.frame(meta, PC1 = M$PCs[, 1], PC2 = M$PCs[, 2])
  p <- ggplot(d, aes(PC1, PC2, colour = Population.assignment)) + geom_point(alpha = .8) +
    labs(x = sprintf("PC1 (%.2f%%)", M$pc_var[1]), y = sprintf("PC2 (%.2f%%)", M$pc_var[2]), colour = "Group") + theme_bw()
  ggsave(fig("Fig1A_PCA.png"), p, width = 7, height = 5, dpi = 300) })

# ---------------- Fig. 3: significant markers by chromosome and variable type ----------------
block("Fig3", {
  d <- hits_all %>% mutate(type = vtype(trait), Mb = POS / 1e6,
         set = factor(set, levels = c("full", "rotundata", "praehensilis"),
                      labels = c("Full panel (Q+K)", "Within D. rotundata", "Within D. praehensilis")))
  p <- ggplot(d, aes(Mb, CHROM, colour = type, shape = type)) + geom_point(size = 2, alpha = .8, position = position_jitter(width = 0, height = .15, seed = 1)) +
    facet_wrap(~set, nrow = 1) + scale_colour_manual(values = c(precipitation = "#2166ac", temperature = "#b2182b", soil = "#1b7837", `climate IC` = "grey30")) +
    labs(x = "Position (Mb)", y = NULL, colour = "Variable type", shape = "Variable type") + theme_bw()
  p3 <- p; ggsave(fig("Fig3_markers_by_chromosome.png"), p, width = 14, height = 6, dpi = 300, bg = "white") })

# ---------------- Fig. 4: A partial RDA, B offset, C-D Köppen ----------------
block("Fig4", {
  need("terra"); out <- readRDS(res("rda_models.rds")); mod <- out$full$penv
  ve <- round(100 * mod$CCA$eig[1:2] / mod$tot.chi, 1)
  si <- as.data.frame(vegan::scores(mod, 1:2, "sites")); si$group <- meta$Population.assignment[match(rownames(si), meta$SampleID)]
  bp <- as.data.frame(vegan::scores(mod, 1:2, "bp")) * 3; bp$var <- rownames(bp)
  bp_lab <- bp[order(-(bp$RDA1^2 + bp$RDA2^2)), ][1:min(10, nrow(bp)), ]   # label the 10 longest arrows
  pA <- ggplot(si, aes(RDA1, RDA2, colour = group)) + geom_point(alpha = .8, size = 1.5) +
    geom_segment(data = bp, aes(0, 0, xend = RDA1, yend = RDA2), inherit.aes = FALSE, colour = "red3", arrow = arrow(length = unit(.12, "cm"))) +
    { if (requireNamespace("ggrepel", quietly = TRUE))
        ggrepel::geom_label_repel(data = bp_lab, aes(RDA1, RDA2, label = var), inherit.aes = FALSE, size = 2.4,
                                  label.size = 0, fill = scales::alpha("white", .8), max.overlaps = Inf, min.segment.length = 0, seed = 1)
      else geom_text(data = bp_lab, aes(RDA1, RDA2, label = var), inherit.aes = FALSE, size = 2.4, check_overlap = TRUE) } +
    labs(tag = "A", x = paste0("RDA1 (", ve[1], "% of total)"), y = paste0("RDA2 (", ve[2], "% of total)"), colour = NULL) + theme_bw()
  r <- terra::rast(res(sprintf("offset_ssp%s_%s.tif", CFG$ssp, CFG$period)))
  od <- as.data.frame(r, xy = TRUE, na.rm = FALSE); od <- od[!is.na(od$extrapolation_masked), ]
  ex <- c(range(od$x), range(od$y))
  pB <- ggplot(od, aes(x, y, fill = offset_mean)) + geom_raster() + scale_fill_viridis_c(name = "Offset", na.value = "grey80") +
    world_borders(xlim = ex[1:2], ylim = ex[3:4], colour = "grey25", linewidth = .3) +
    geom_point(data = meta, aes(Longitude, Latitude), inherit.aes = FALSE, shape = 21, fill = "white", size = 1) +
    coord_quickmap(xlim = ex[1:2], ylim = ex[3:4], expand = FALSE) + labs(tag = "B", x = NULL, y = NULL) + theme_bw()
  kg <- c("Af","Am","Aw","BWh","BWk","BSh","BSk","Csa","Csb","Csc","Cwa","Cwb","Cwc","Cfa","Cfb","Cfc",
          "Dsa","Dsb","Dsc","Dsd","Dwa","Dwb","Dwc","Dwd","Dfa","Dfb","Dfc","Dfd","ET","EF")
  kcol <- c(Af = "#0000FF", Am = "#0078FF", Aw = "#46AAFA", BWh = "#FF0000", BWk = "#FF9696", BSh = "#F5A500", BSk = "#FFDC64",
            Cwa = "#96FF96", Cwb = "#64C864", Cfa = "#C8FF50", Cfb = "#64FF50", Csa = "#FFFF00", Csb = "#C8C800")
  kplot <- function(f, ttl) { k <- terra::aggregate(terra::crop(terra::rast(f), terra::ext(ex)), 5, fun = "modal")
    kd <- as.data.frame(k, xy = TRUE); names(kd)[3] <- "code"; kd$class <- factor(kg[kd$code])
    ggplot(kd, aes(x, y, fill = class)) + geom_raster() + scale_fill_manual(values = kcol, na.value = "grey85", drop = FALSE) +
      world_borders(xlim = ex[1:2], ylim = ex[3:4], colour = "grey25", linewidth = .3) +
      geom_point(data = meta, aes(Longitude, Latitude), inherit.aes = FALSE, shape = 21, fill = "white", size = 1) +
      coord_quickmap(xlim = ex[1:2], ylim = ex[3:4], expand = FALSE) + labs(tag = ttl, x = NULL, y = NULL, fill = "Köppen") + theme_bw() }
  pC <- kplot(CFG$koppen_hist, "C"); pD <- kplot(CFG$koppen_fut, "D")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    fig4 <- patchwork::wrap_plots(pA, pB, pC, pD, ncol = 2); ggsave(fig("Fig4_composite.png"), fig4, width = 15, height = 11, dpi = 300, bg = "white")
  } else stop("patchwork is required for the Figure 4 composite") })

# ---------------- Fig. S4: FST outliers by chromosome ----------------
block("FigS4", {
  d <- read.csv(tab("FigS4_outliers_by_chromosome.csv")); d$outlier <- factor(d$outlier, c("95th", "99th", "99.5th"))
  p <- ggplot(d, aes(CHROM, n, fill = outlier)) + geom_col() + scale_fill_manual(values = c("#ffffb2", "#fd8d3c", "#bd0026")) +
    coord_flip() + labs(x = NULL, y = "Outlier SNPs (global WC FST)", fill = "Percentile") + theme_bw()
  ggsave(fig("FigS4_FST_outliers_by_chromosome.png"), p, width = 6, height = 6, dpi = 300) })

# ---------------- Fig. S5 / S6: environmental correlations ----------------
block("FigS5", {
  R <- cor(M$E); o <- hclust(as.dist(1 - abs(R)))$order
  d <- as.data.frame(as.table(R[o, o]))
  p <- ggplot(d, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient2(low = "red", mid = "white", high = "blue", limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, fill = "r") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, size = 6), axis.text.y = element_text(size = 6))
  ggsave(fig("FigS5_env_correlation_matrix.png"), p, width = 8, height = 7, dpi = 300) })
block("FigS6", {
  need("igraph"); R <- cor(M$E); A <- (abs(R) >= .3) * 1; diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected"); set.seed(CFG$seed); lay <- igraph::layout_with_fr(g)
  el <- igraph::as_edgelist(g); ed <- data.frame(x = lay[match(el[, 1], colnames(R)), 1], y = lay[match(el[, 1], colnames(R)), 2],
        xend = lay[match(el[, 2], colnames(R)), 1], yend = lay[match(el[, 2], colnames(R)), 2], r = R[el])
  nd <- data.frame(x = lay[, 1], y = lay[, 2], var = colnames(R), type = vtype(colnames(R)))
  p <- ggplot() + geom_segment(data = ed, aes(x, y, xend = xend, yend = yend, colour = r), linewidth = .4) +
    scale_colour_gradient2(low = "red", mid = "grey90", high = "darkgreen", limits = c(-1, 1)) +
    geom_point(data = nd, aes(x, y, fill = type), shape = 21, size = 5) + geom_text(data = nd, aes(x, y, label = var), size = 2.2, vjust = -1.4) +
    scale_fill_manual(values = c(precipitation = "#4393c3", temperature = "#d6604d", soil = "#5aae61")) + theme_void() + labs(colour = "r", fill = NULL)
  ggsave(fig("FigS6_env_correlation_network.png"), p, width = 9, height = 8, dpi = 300, bg = "white") })

# ---------------- Fig. S7: EGWAS summary (full panel) ----------------
block("FigS7", {
  a <- hits %>% count(trait) %>% ggplot(aes(reorder(trait, n), n)) + geom_col() + coord_flip() + labs(x = NULL, y = "Significant markers", tag = "A") + theme_bw()
  b <- hits %>% distinct(SNP, CHROM) %>% count(CHROM) %>% ggplot(aes(CHROM, n)) + geom_col() + coord_flip() + labs(x = NULL, y = "Markers", tag = "B") + theme_bw()
  c <- ggplot(hits, aes(-log10(p))) + geom_histogram(bins = 30) + geom_vline(xintercept = median(-log10(hits$p)), colour = "red", linetype = 2) + labs(tag = "C") + theme_bw()
  d <- ggplot(hits, aes(reorder(trait, -log10(p), median), 100 * PVE)) + geom_boxplot() + coord_flip() + labs(x = NULL, y = "PVE (%)", tag = "D") + theme_bw()
  ggsave(fig("FigS7_EGWAS_summary.png"), patchwork::wrap_plots(a, b, c, d, ncol = 2), width = 12, height = 10, dpi = 300, bg = "white") })

# ---------------- Fig. S8 / S9 / S11 / S13: pleiotropic markers ----------------
P <- if (file.exists(res("egwas_pvalues_full.rds"))) readRDS(res("egwas_pvalues_full.rds")) else NULL
top_multi <- if (!is.null(hits)) hits %>% count(SNP, sort = TRUE) %>% filter(n > 1) else NULL
block("FigS8", {
  sn <- head(top_multi$SNP, 30); d <- hits %>% filter(SNP %in% sn)
  p <- ggplot(d, aes(trait, factor(SNP, levels = rev(sn)), fill = -log10(p))) + geom_tile() + scale_fill_viridis_c() +
    labs(x = NULL, y = NULL, fill = "-log10 p") + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
  ggsave(fig("FigS8_pleiotropic_heatmap.png"), p, width = 9, height = 8, dpi = 300) })
block("FigS9", {
  sn <- head(top_multi$SNP, 15); d <- hits %>% filter(SNP %in% sn); vars <- unique(d$trait)
  ang_v <- seq(0, 2 * pi, length = length(vars) + 1)[-1]; ang_s <- seq(0, 2 * pi, length = length(sn) + 1)[-1]
  V <- data.frame(trait = vars, x = 2 * cos(ang_v), y = 2 * sin(ang_v)); S_ <- data.frame(SNP = sn, xs = cos(ang_s), ys = sin(ang_s))
  e <- d %>% left_join(V, by = "trait") %>% left_join(S_, by = "SNP")
  p <- ggplot() + geom_segment(data = e, aes(xs, ys, xend = x, yend = y, colour = SNP, linewidth = -log10(p)), alpha = .6) +
    geom_point(data = S_, aes(xs, ys, colour = SNP), size = 3) + geom_label(data = V, aes(x, y, label = trait), size = 2.5) +
    geom_text(data = S_, aes(xs * 1.25, ys * 1.25, label = sub("^BDMI0*", "", SNP), colour = SNP), size = 2.6, fontface = "bold") +
    scale_linewidth(range = c(.2, 1.5)) + coord_equal(clip = "off") + theme_void() + guides(colour = "none")
  ggsave(fig("FigS9_pleiotropy_network.png"), p, width = 8, height = 8, dpi = 300, bg = "white") })
block("FigS11", {
  sn <- head(top_multi$SNP, 15)
  d <- hits %>% filter(SNP %in% sn) %>% group_by(SNP) %>% slice_min(p, n = 3) %>% ungroup()
  long <- bind_rows(lapply(seq_len(nrow(d)), function(i) data.frame(SNP = d$SNP[i], trait = d$trait[i],
            genotype = factor(M$G[, d$SNP[i]], 0:2, c("Ref/Ref", "Het", "Alt/Alt")), value = Y[rownames(M$G), d$trait[i]])))
  long <- long[!is.na(long$genotype), ]
  p <- ggplot(long, aes(genotype, value, fill = genotype)) + geom_boxplot(outlier.size = .5) +
    facet_wrap(~ SNP + trait, scales = "free", ncol = 5, labeller = labeller(.multi_line = TRUE)) + scale_fill_manual(values = c("Ref/Ref" = "#4575b4", "Het" = "grey70", "Alt/Alt" = "#d73027")) +
    scale_x_discrete(drop = TRUE) + labs(x = NULL, y = "Environmental value") + theme_bw() + theme(strip.text = element_text(size = 8), axis.text = element_text(size = 7), legend.position = "none")
  ggsave(fig("FigS11_env_by_genotype_top_markers.png"), p, width = 13, height = 2.6 * ceiling(nrow(d) / 5) + 1, dpi = 300, limitsize = FALSE) })
block("FigS13", {
  sn <- intersect(head(top_multi$SNP, 5), colnames(M$G))
  if (!length(sn)) stop("no multi-variable markers found among genotyped SNPs")
  d <- bind_rows(lapply(sn, function(s) data.frame(meta[, c("Longitude", "Latitude", "Population.assignment")], SNP = s,
          genotype = factor(M$G[, s], 0:2, c("Ref/Ref", "Het", "Alt/Alt")))))
  p <- ggplot(d, aes(Longitude, Latitude)) + world_borders(xlim = range(d$Longitude) + c(-1, 1), ylim = range(d$Latitude) + c(-1, 1), colour = "grey60") +
    geom_point(aes(colour = genotype, shape = genotype), position = position_jitter(.08, .08), size = 1.8) +
    scale_colour_manual(values = c("#4575b4", "grey50", "#d73027")) + facet_wrap(~SNP) +
    coord_quickmap(xlim = range(d$Longitude) + c(-.5, .5), ylim = range(d$Latitude) + c(-.5, .5)) +
    labs(x = "Longitude", y = "Latitude", colour = "Genotype", shape = "Genotype") + theme_bw()
  ggsave(fig("FigS13_allele_maps_multivariable_markers.png"), p, width = 12, height = 3.6 * ceiling(length(sn) / 3), dpi = 300, bg = "white")
  write.csv(data.frame(SNP = sn), tab("FigS13_marker_ids.csv"), row.names = FALSE) })   # use these IDs in the caption


# ---------------- Fig. S15: within-group GEAV heatmap (redrawn from Table S6, no refitting) ----------------
block("FigS15", {
  out <- read.csv(tab("TableS6_GEAV_within_group.csv"), check.names = FALSE)
  cp <- if (file.exists(tab("candidate_parents_within_group.csv"))) read.csv(tab("candidate_parents_within_group.csv"))$SampleID else character(0)
  pS15 <- plot_geav_heatmap(out, M, cand_ids = cp)
  sc <- attr(pS15, "subclusters"); write.csv(sc, tab("TableS6b_GEAV_subclusters.csv"), row.names = FALSE)
  print(sc %>% filter(!is.na(sep_ratio)) %>% group_by(group) %>% summarise(n = n(), sep_ratio = round(first(sep_ratio), 2),
        sizes = paste(table(cluster), collapse = "/"), .groups = "drop"))
  ggsave(fig("FigS15_GEAV_within_group.png"), pS15, width = 10, height = 16, dpi = 300, bg = "white")
})

# ---------------- Fig. S16: current vs future climate at sampling sites ----------------
block("FigS16", {
  cs <- readRDS(res("climate_at_sites.rds")); now <- cs$now; fut <- cs$fut; colnames(fut) <- colnames(now)
  pc <- prcomp(rbind(now, fut), scale. = TRUE); ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1); n <- nrow(now)
  d <- data.frame(x0 = pc$x[1:n, 1], y0 = pc$x[1:n, 2], x1 = pc$x[n + 1:n, 1], y1 = pc$x[n + 1:n, 2], group = cs$group)
  p <- ggplot(d) + geom_segment(aes(x0, y0, xend = x1, yend = y1, colour = group), arrow = arrow(length = unit(.08, "cm")), alpha = .6) +
    geom_point(aes(x0, y0, colour = group), size = 1) +
    labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"), colour = NULL) + theme_bw()
  ggsave(fig("FigS16_climate_shift_at_sites.png"), p, width = 8, height = 6, dpi = 300)
  write.csv(data.frame(PC = 1:3, var_pct = ve[1:3]), tab("FigS16_PC_variance.csv"), row.names = FALSE) })

# ---------------- Composite supplementary figures (multi-part) ----------------
# Panels are read from results/<set>/panels/ (saved by scripts 02, 04 and 07) and assembled here.
need("patchwork")
block("FigS1 composite", {          # built from the master object (no dependency on script 02 panels)
  chr_n <- dplyr::count(M$snp, CHROM); chr_n$chr <- factor(chrom_number(chr_n$CHROM), levels = sort(unique(chrom_number(chr_n$CHROM))))
  maf_g <- bind_rows(lapply(split(meta$SampleID, meta$Population.assignment), function(ids) {
    p <- colMeans(M$G[ids, , drop = FALSE], na.rm = TRUE) / 2; data.frame(maf = pmin(p, 1 - p)) }), .id = "group")
  a <- ggplot(chr_n, aes(chr, n)) + geom_col(fill = "grey40") + coord_flip() + labs(x = "Pseudomolecule", y = "SNPs") + theme_bw()
  b <- ggplot(maf_g, aes(maf)) + geom_histogram(bins = 25, fill = "grey40") + facet_wrap(~group, scales = "free_y") +
    labs(x = "Minor allele frequency (within group)", y = "SNPs") + theme_bw() + theme(strip.text = element_text(size = 8))
  f <- (a | b) + patchwork::plot_layout(widths = c(1, 2.2)) + patchwork::plot_annotation(tag_levels = "A")
  ggsave(fig("FigS1_composite.png"), f, width = 13, height = 7, dpi = 300, bg = "white") })
# Figures S2 and S3 from the TESS3 run (Yam_spatial_population_structure_revised.R -> results/strict/tess3_K2_7.rds)
t3f <- file.path("results", "strict", "tess3_K2_7.rds")
block("FigS2 composite", {
  if (!file.exists(t3f)) stop("results/strict/tess3_K2_7.rds not found - run Yam_spatial_population_structure_revised.R")
  t3 <- readRDS(t3f); grp <- ifelse(is.na(t3$group), "Unassigned", t3$group)
  ordK <- order(grp, apply(t3$Q$K7, 1, which.max), -apply(t3$Q$K7, 1, max))     # same accession order in every panel
  bars <- lapply(names(t3$Q), function(k) {
    Q <- t3$Q[[k]]; K <- ncol(Q)
    d <- data.frame(Sample = factor(t3$sample, levels = t3$sample[ordK]), group = grp, Q)
    d <- tidyr::pivot_longer(d, -c(Sample, group), names_to = "cluster", values_to = "q")
    d$cluster <- factor(d$cluster, levels = colnames(data.frame(Q)))
    ggplot(d, aes(Sample, q, fill = cluster)) + geom_col(width = 1, colour = NA) +
      facet_grid(~group, scales = "free_x", space = "free_x", labeller = label_wrap_gen(width = 14)) +
      scale_fill_manual(values = unname(t3$bar_palettes[[k]][1:K]), guide = "none") +
      scale_y_continuous(expand = c(0, 0), breaks = c(0, .5, 1)) +
      labs(x = NULL, y = paste0("K = ", K)) + theme_classic() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(1, "pt"),
            strip.text = if (k == "K2") element_text(size = 6.5, angle = 0) else element_blank(),
            strip.background = element_blank())
  })
  f <- patchwork::wrap_plots(bars, ncol = 1) + patchwork::plot_annotation(tag_levels = "A")
  ggsave(fig("FigS2_composite.png"), f, width = 13, height = 10, dpi = 300, bg = "white") })
block("FigS3 composite", {
  if (!file.exists(t3f)) stop("results/strict/tess3_K2_7.rds not found - run Yam_spatial_population_structure_revised.R")
  if (!requireNamespace("tess3r", quietly = TRUE)) stop("install tess3r to draw Figure S3")
  t3 <- readRDS(t3f); xy <- t3$coords
  xl <- range(xy[, 1]) + c(-1.5, 1.5); yl <- range(xy[, 2]) + c(-1.5, 1.5)
  maps <- lapply(names(t3$Q), function(k)
    tess3r::ggtess3Q(t3$Q[[k]], xy, map.polygon = t3$map.polygon, col.palette = t3$map_palettes[[k]]) +
      geom_path(data = t3$map.polygon, aes(x = long, y = lat, group = group), linewidth = .25, colour = "grey30") +
      geom_point(data = data.frame(Longitude = xy[, 1], Latitude = xy[, 2]), aes(Longitude, Latitude), size = .4) +
      coord_equal(xlim = xl, ylim = yl, expand = FALSE) +
      labs(subtitle = paste0("K = ", ncol(t3$Q[[k]])), x = "Longitude", y = "Latitude") + theme_bw() +
      theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9)))
  f <- patchwork::wrap_plots(maps, ncol = 2) + patchwork::plot_annotation(tag_levels = "A")
  ggsave(fig("FigS3_composite.png"), f, width = 12, height = 13, dpi = 300, bg = "white") })
block("FigS10 composite", {
  pf <- res("egwas_pvalues_full.rds")
  rows <- lapply(paste0("IC", 1:3), function(tr) {
    pp <- if (file.exists(pf)) egwas_manhattan_qq(readRDS(pf)[[tr]], M$snp) else
      list(manhattan = load_panel(paste0("FigS10_", tr, "_manhattan"), "07_EGWAS_revised.R"),
           qq = load_panel(paste0("FigS10_", tr, "_qq"), "07_EGWAS_revised.R"))
    if (!is.null(pp$lambda)) msg(tr, ": lambda =", round(pp$lambda, 2), "| SNPs at FDR", CFG$fdr, "=", pp$n_sig)
    patchwork::wrap_elements(full = pp$manhattan + pp$qq + patchwork::plot_layout(widths = c(3, 1)))
  })
  f <- patchwork::wrap_plots(rows, ncol = 1) + patchwork::plot_annotation(tag_levels = "A")
  ggsave(fig("FigS10_composite.png"), f, width = 14, height = 11, dpi = 300, bg = "white") })
block("FigS13 composite", {
  sn <- intersect(head(top_multi$SNP, 5), colnames(M$G))
  if (!length(sn)) stop("no multi-variable markers found among genotyped SNPs")
  xl <- range(meta$Longitude, na.rm = TRUE); yl <- range(meta$Latitude, na.rm = TRUE)
  maps13 <- lapply(sn, function(s) {
    d <- data.frame(meta[, c("Longitude", "Latitude")], genotype = factor(M$G[, s], 0:2, c("Ref/Ref", "Het", "Alt/Alt")))
    ggplot(d, aes(Longitude, Latitude)) + world_borders(xlim = xl + c(-1, 1), ylim = yl + c(-1, 1), colour = "grey60") +
      geom_point(aes(colour = genotype, shape = genotype), position = position_jitter(.08, .08, seed = 1), size = 1.8) +
      scale_colour_manual(values = c("Ref/Ref" = "#4575b4", "Het" = "grey50", "Alt/Alt" = "#d73027"), drop = FALSE) +
      scale_shape_manual(values = c("Ref/Ref" = 16, "Het" = 17, "Alt/Alt" = 15), drop = FALSE) +
      coord_quickmap(xlim = xl + c(-.5, .5), ylim = yl + c(-.5, .5)) +
      labs(subtitle = s, x = "Longitude", y = "Latitude", colour = "Genotype", shape = "Genotype") + theme_bw() +
      theme(plot.subtitle = element_text(size = 8))
  })
  f <- patchwork::wrap_plots(maps13, nrow = 1) + patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")
  ggsave(fig("FigS13_composite.png"), f, width = 4.2 * length(sn) + 1.5, height = 3.6, dpi = 300, bg = "white") })
block("FigS14 composite", {
  acc <- read.csv(tab("FigS14_prediction_accuracy.csv"))
  schemes <- c("random 10-fold", "leave-one-group-out", "spatial blocks")
  rng <- range(acc$r_pooled, na.rm = TRUE)
  pl <- lapply(schemes, function(sc) ggplot(acc[acc$scheme == sc, ], aes(trait, r_pooled, colour = model)) +
    geom_point(position = position_dodge(.6), size = 1.3) + geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(limits = rng) + labs(x = NULL, y = "Prediction accuracy (r)", colour = "Model", subtitle = sc) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6, vjust = .5, hjust = 1), plot.subtitle = element_text(size = 9)))
  f <- patchwork::wrap_plots(pl, ncol = 1) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(tag_levels = "A")
  ggsave(fig("FigS14_composite.png"), f, width = 12, height = 11, dpi = 300, bg = "white") })
block("FigS17 composite", {
  a <- load_panel("FigS17A", "04_environment_IBD_IBE.R"); b <- load_panel("FigS17B", "04_environment_IBD_IBE.R")
  f <- (a | b) + patchwork::plot_annotation(tag_levels = "A")
  ggsave(fig("FigS17_composite.png"), f, width = 15, height = 6, dpi = 300, bg = "white") })

# ---------------- Main figures for submission: separate high-resolution files ----------------
# Writes Figure1.tiff ... Figure4.tiff (600 dpi, LZW) and matching PDFs to figures/<set>/submission/.
# Figure 1B-C (TESS3 maps) come from Yam_spatial_population_structure.R: if
# figures/<set>/Fig1BC_tess3.png exists it is added as a panel; otherwise Figure 1 is panel A only.
block("Main figures", {
  need("patchwork"); sub_dir <- file.path(FIG, "submission"); dir.create(sub_dir, showWarnings = FALSE)
  save_main <- function(p, name, w, h) {
    ggsave(file.path(sub_dir, paste0(name, ".tiff")), p, width = w, height = h, dpi = 600, compression = "lzw", bg = "white")
    ggsave(file.path(sub_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white") }
  # Figure 1: A PCA (this script); B TESS3 pie map and C interpolated ancestry (K = 7) from
  # Yam_spatial_population_structure_revised.R (results/strict/Fig1BC_panels.rds)
  d <- data.frame(meta, PC1 = M$PCs[, 1], PC2 = M$PCs[, 2])
  f1a <- ggplot(d, aes(PC1, PC2, colour = Population.assignment)) + geom_point(alpha = .8, size = 1.6) +
    labs(x = sprintf("PC1 (%.1f%%)", M$pc_var[1]), y = sprintf("PC2 (%.1f%%)", M$pc_var[2]), colour = "Genetic group") +
    theme_bw() + theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8))
  rds <- file.path("results", "strict", "Fig1BC_panels.rds")
  if (file.exists(rds)) {
    for (pk in c("tess3r", "scatterpie")) if (!requireNamespace(pk, quietly = TRUE)) stop("Install ", pk, " to draw Figure 1B-C")
    suppressPackageStartupMessages({ library(scatterpie); library(tess3r) })
    bc <- readRDS(rds)
    bottom <- (bc$B + bc$C) + patchwork::plot_layout(guides = "collect")     # shared ancestry legend at right
    f1 <- patchwork::wrap_elements(full = f1a + theme(aspect.ratio = .75)) / bottom +
      patchwork::plot_layout(heights = c(1.1, 1)) + patchwork::plot_annotation(tag_levels = "A")
    w1 <- 11; h1 <- 8.5
  } else if (file.exists(fig("Fig1BC_tess3.png")) && requireNamespace("png", quietly = TRUE)) {
    f1 <- f1a / patchwork::wrap_elements(full = grid::rasterGrob(png::readPNG(fig("Fig1BC_tess3.png")))) +
      patchwork::plot_annotation(tag_levels = "A"); w1 <- 10; h1 <- 9
  } else { msg("Figure 1: no TESS3 panels found - run Yam_spatial_population_structure_revised.R first; panel A only"); f1 <- f1a; w1 <- 7; h1 <- 5 }
  ggsave(fig("Fig1_composite.png"), f1, width = w1, height = h1, dpi = 300, bg = "white")
  save_main(f1, "Figure1", w1, h1)
  # Figure 2 composite (A FST by pseudomolecule, B rarefied and C raw private alleles, D FROH) from saved tables
  f2 <- fig2_composite(read.csv(res("fst_per_locus_WC.csv")), read.csv(tab("Fig2A_private_alleles_by_group.csv")),
                       read.csv(tab("Fig2B_froh_by_individual.csv")))
  ggsave(fig("Fig2_composite.png"), f2, width = 11, height = 7.5, dpi = 300, bg = "white")
  save_main(f2, "Figure2", 11, 7.5)
  # Figure 3 and Figure 4: the PNG versions drawn above are re-saved at 600 dpi if their plot objects exist
  for (fx in list(list("p3", "Figure3", 10, 5.5), list("fig4", "Figure4", 10, 7.5)))
    if (exists(fx[[1]])) tryCatch(save_main(get(fx[[1]]), fx[[2]], fx[[3]], fx[[4]]),
                                  error = function(e) msg("SKIPPED", fx[[2]], "-", conditionMessage(e)))
  msg("Submission figures written to", sub_dir)
})
