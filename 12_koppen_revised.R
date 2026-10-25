################################################################################
# 12_koppen_revised.R                                      Analysis Plan: A8 (Köppen)
# Short replacement for the extraction/summary part of
# koppen_geiger_climate_change_analysis.R. The original reads the 2071_2099
# raster but labels every plot and table "2041-2070" -> labels corrected here.
# The mapping figures of the original can be kept; change their subtitle text to
# "1991-2020 to 2071-2099 (SSP5-8.5)".
# Output: Köppen class per site (now / future), change table, link to offsets.
################################################################################
source("00_config_and_helpers.R"); need("terra"); M <- load_master(); meta <- M$meta
kg <- c("Af","Am","Aw","BWh","BWk","BSh","BSk","Csa","Csb","Csc","Cwa","Cwb","Cwc","Cfa","Cfb","Cfc",
        "Dsa","Dsb","Dsc","Dsd","Dwa","Dwb","Dwc","Dwd","Dfa","Dfb","Dfc","Dfd","ET","EF")   # Beck et al. 2023 codes 1-30
pts <- terra::vect(meta, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
meta$KG_1991_2020 <- kg[terra::extract(terra::rast(CFG$koppen_hist), pts)[, 2]]
meta$KG_2071_2099 <- kg[terra::extract(terra::rast(CFG$koppen_fut), pts)[, 2]]
sites <- meta %>% distinct(SiteID, .keep_all = TRUE)
tr <- sites %>% count(KG_1991_2020, KG_2071_2099, name = "n_sites") %>% mutate(changed = KG_1991_2020 != KG_2071_2099)
write.csv(tr, tab("Koppen_transitions_sites.csv"), row.names = FALSE); print(tr)
msg("Sites changing class:", sum(sites$KG_1991_2020 != sites$KG_2071_2099), "of", nrow(sites))
f_off <- tab("TableS6_offset_per_accession.csv")
if (file.exists(f_off)) {
  off <- read.csv(f_off) %>% left_join(meta[, c("SampleID", "KG_1991_2020", "KG_2071_2099")], by = "SampleID") %>%
    mutate(changed = KG_1991_2020 != KG_2071_2099)
  print(off %>% group_by(changed) %>% summarise(n = n(), median_offset = median(offset_mean)))
  if (length(unique(off$changed)) == 2) print(wilcox.test(offset_mean ~ changed, data = off))
  write.csv(off, tab("TableS6_offset_with_Koppen.csv"), row.names = FALSE)
}
