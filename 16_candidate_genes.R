################################################################################
# 16_candidate_genes.R                                     Discussion (candidate genes)
# No original script uploaded (genes were looked up manually in the Yam Genome Hub).
# Lists genes within +/- window of the revised candidate markers, prioritising markers
# that (i) are significant within species, (ii) fall in top WZA windows, or (iii) are
# shared by EGWAS and LFMM. Optional BLASTP of their proteins against TAIR10.
# Inputs: CFG$gff (TDr96 GFF3; chromosome names must match the VCF), optional
#         CFG$proteins_fa and CFG$tair_db with blastp on the PATH.
# Outputs: TableS11_candidate_markers_prioritised.csv, TableS12_candidate_genes.csv, blast_TAIR10.tsv
################################################################################
source("00_config_and_helpers.R"); M <- load_master()
hits <- read.csv(tab("TableS18_EGWAS_significant_with_effects.csv"))
wza  <- if (file.exists(tab("TableS6_WZA_top_windows.csv"))) read.csv(tab("TableS6_WZA_top_windows.csv")) else NULL
lfmm <- if (file.exists(tab("TableS8_LFMM_significant.csv"))) read.csv(tab("TableS8_LFMM_significant.csv")) else NULL
ld   <- if (file.exists(tab("LD_decay_summary.csv"))) read.csv(tab("LD_decay_summary.csv")) else NULL
win  <- max(CFG$gene_window, if (!is.null(ld)) ld$half_decay_bp[ld$set == "rotundata"] else 0, na.rm = TRUE)

mk <- hits %>% group_by(SNP, CHROM, POS) %>%
  summarise(n_vars = n_distinct(trait[set == "full"]), within_species = any(set != "full"),
            sets = paste(sort(unique(set)), collapse = ";"), traits = paste(sort(unique(trait)), collapse = ";"),
            max_PVE = max(PVE), min_q = min(q), .groups = "drop")
mk$in_WZA_top <- if (!is.null(wza)) mapply(function(ch, p) any(wza$chr == ch & p >= wza$start & p < wza$start + CFG$wza_window), mk$CHROM, mk$POS) else NA
mk$LFMM_shared <- if (!is.null(lfmm)) mk$SNP %in% lfmm$SNP[lfmm$shared_with_EGWAS] else NA
mk$priority <- mk$within_species * 2 + coalesce(mk$in_WZA_top, FALSE) + coalesce(mk$LFMM_shared, FALSE) + (mk$n_vars > 1)
mk <- arrange(mk, desc(priority), min_q)
write.csv(mk, tab("TableS11_candidate_markers_prioritised.csv"), row.names = FALSE)

if (!file.exists(CFG$gff)) stop("GFF not found: set CFG$gff. Marker list written to TableS11_candidate_markers_prioritised.csv")
gff <- fread(cmd = paste("grep -v '^#'", shQuote(CFG$gff)), sep = "\t", header = FALSE, quote = "",
             col.names = c("seqid", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
genes <- gff[type %in% c("gene", "mRNA")]
shared <- intersect(unique(genes$seqid), unique(M$snp$CHROM))
msg("GFF:", CFG$gff, "|", nrow(genes), "gene/mRNA features on", length(unique(genes$seqid)), "sequences;",
    length(shared), "sequence names match the VCF")
if (!length(shared)) stop("No GFF sequence names match the VCF chromosome names (VCF e.g. ", M$snp$CHROM[1],
                          "; GFF e.g. ", genes$seqid[1], "). Rename the GFF seqids to the VCF names first.")
genes[, ID := sub(".*ID=([^;]+).*", "\\1", attr)]
genes[, product := ifelse(grepl("product=", attr), sub(".*product=([^;]+).*", "\\1", attr), NA)]
genes[, name := ifelse(grepl("Name=", attr), sub(".*Name=([^;]+).*", "\\1", attr), NA)]
top <- mk %>% filter(priority >= 2)
if (!nrow(top)) top <- head(mk, 20)
cg <- rbindlist(lapply(seq_len(nrow(top)), function(i) {
  g <- genes[seqid == top$CHROM[i] & end >= top$POS[i] - win & start <= top$POS[i] + win]
  if (!nrow(g)) return(NULL)
  g[, `:=`(SNP = top$SNP[i], dist_bp = pmax(0, pmax(start, top$POS[i]) - pmin(end, top$POS[i])), traits = top$traits[i], sets = top$sets[i])]
  g[, .(SNP, traits, sets, seqid, start, end, type, ID, name, product, dist_bp)] }))
write.csv(cg, tab("TableS12_candidate_genes.csv"), row.names = FALSE)
msg("Markers searched:", nrow(top), "; window +/-", win, "bp; gene/mRNA features found:", nrow(cg))
if (nrow(cg) && file.exists(CFG$proteins_fa) && nzchar(Sys.which("blastp"))) {
  writeLines(unique(cg$ID), res("candidate_ids.txt"))
  system(sprintf("awk 'BEGIN{while((getline l < \"%s\")>0) k[l]=1} /^>/{split(substr($1,2),a,\" \"); p=(a[1] in k)} p' %s > %s",
                 res("candidate_ids.txt"), CFG$proteins_fa, res("candidate_proteins.fa")))
  system(sprintf("blastp -query %s -db %s -evalue 1e-10 -max_target_seqs 3 -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' > %s",
                 res("candidate_proteins.fa"), CFG$tair_db, tab("blast_TAIR10.tsv")))
  msg("BLAST results:", tab("blast_TAIR10.tsv"), "(protein FASTA headers must use the same IDs as the GFF)")
}
