################################################################################
# 03_relaxed_missing_dataset.R                             Analysis Plan: A2 (part 1)
#
# Creates the relaxed SNP set (<= 20% missing, MAF >= 0.05, biallelic SNPs),
# imputes it with Beagle, then you rebuild the master object:
#     YAM_GENO_SET=relaxed Rscript 01_build_master_data.R
# and rerun 05 and 07 with YAM_GENO_SET=relaxed. Script 19 compares the two sets.
#
# Requires bcftools, tabix/bgzip (htslib) and Java.
# Software paths can be set here or with environment variables
#   YAM_JAVA, YAM_BEAGLE_JAR, YAM_BCFTOOLS, YAM_REF_FAI
#
# Fixes in this version
#   - Full paths to java / bcftools (RStudio on macOS often does not inherit the
#     shell PATH, so system("bcftools ...") can fail even when it works in Terminal).
#   - Adds ##contig lines to the VCF header. The original VCF has none, which
#     causes "[W::vcf_parse] Contig '...' is not defined in the header" and can
#     make indexing or Beagle fail. Lengths come from the reference .fai if given
#     (preferred), otherwise from the largest observed position per contig.
#   - Checks the input is bgzip-compressed (plain gzip cannot be indexed).
#   - Indexes the Beagle output and checks that SNP and sample counts are sane.
################################################################################
source("00_config_and_helpers.R")

# ---- software paths ------------------------------------------------------------
java_bin   <- Sys.getenv("YAM_JAVA",       "/opt/homebrew/opt/openjdk/bin/java")
beagle_jar <- path.expand(Sys.getenv("YAM_BEAGLE_JAR", "~/software/beagle/beagle.29Oct24.c8e.jar"))
bcftools   <- Sys.getenv("YAM_BCFTOOLS",   Sys.which("bcftools"))
if (!nzchar(bcftools)) bcftools <- "/opt/homebrew/bin/bcftools"
bgzip      <- file.path(dirname(bcftools), "bgzip"); if (!file.exists(bgzip)) bgzip <- Sys.which("bgzip")
ref_fai    <- Sys.getenv("YAM_REF_FAI", "")   # e.g. "TDr96_F1_v2_PseudoChromosome.fasta.fai" (optional)
threads    <- 8; java_mem <- "-Xmx16g"

for (f in c(java_bin, beagle_jar, bcftools, CFG$vcf_full))
  if (!file.exists(f)) stop("Not found: ", f, "  (set the path at the top of this script)")
run <- function(cmd) { msg(cmd); rc <- system(cmd); if (rc != 0) stop("Command failed (exit ", rc, "): ", cmd) }
sq <- shQuote

inp        <- CFG$vcf_full
fixed      <- "Yam_full_with_contigs.vcf.gz"
tmp        <- "Yam_relaxed_unimputed.vcf.gz"
out_prefix <- sub("\\.vcf\\.gz$", "", CFG$vcf_relaxed)

# ---- 0a. input must be BGZF, not plain gzip ---------------------------------------
fmt <- system(paste(sq(file.path(dirname(bcftools), "htsfile")), sq(inp)), intern = TRUE)
if (length(fmt) && grepl("gzip-compressed", fmt) && !grepl("BGZF", fmt)) {
  msg("Input is plain gzip; recompressing with bgzip")
  run(sprintf("gzip -dc %s | %s -@ %d -c > %s", sq(inp), sq(bgzip), threads, sq("Yam_full_bgzf.vcf.gz")))
  inp <- "Yam_full_bgzf.vcf.gz"
}

# ---- 0b. add ##contig lines to the header if missing -------------------------------
n_contig_hdr <- as.integer(system(sprintf("%s view -h %s | grep -c '^##contig' || true", sq(bcftools), sq(inp)), intern = TRUE))
if (is.na(n_contig_hdr) || n_contig_hdr == 0) {
  msg("VCF header has no ##contig lines - adding them")
  if (nzchar(ref_fai) && file.exists(ref_fai)) {
    run(sprintf("%s reheader --fai %s -o %s %s", sq(bcftools), sq(ref_fai), sq(fixed), sq(inp)))
  } else {
    # contig names in file order and maximum position (stand-in for length)
    cp <- data.table::fread(cmd = sprintf("%s query -f '%%CHROM\\t%%POS\\n' %s 2>/dev/null", sq(bcftools), sq(inp)),
                            header = FALSE, col.names = c("chr", "pos"))
    cl <- cp[, .(len = max(pos)), by = chr]
    writeLines(sprintf("##contig=<ID=%s,length=%d>", cl$chr, as.integer(cl$len)), "contig_header_lines.txt")
    msg("Contigs found:", nrow(cl), "(lengths = max observed position; supply YAM_REF_FAI for true lengths)")
    run(sprintf("%s annotate -h contig_header_lines.txt %s -Oz -o %s", sq(bcftools), sq(inp), sq(fixed)))
  }
  run(sprintf("%s index -f %s", sq(bcftools), sq(fixed)))
  inp <- fixed
} else msg("Header already has", n_contig_hdr, "##contig lines")

# ---- 1. filter: biallelic SNPs, <= relaxed missingness, MAF >= threshold -----------
run(sprintf("%s view --threads %d -m2 -M2 -v snps -i 'F_MISSING<=%s && MAF>=%s' %s -Oz -o %s",
            sq(bcftools), threads, CFG$relaxed_missing, CFG$maf, sq(inp), sq(tmp)))
run(sprintf("%s index -f %s", sq(bcftools), sq(tmp)))
n_tmp <- trimws(system(sprintf("%s index -n %s", sq(bcftools), sq(tmp)), intern = TRUE))
msg("SNPs after relaxed filter:", n_tmp)

# ---- 2. impute with Beagle ----------------------------------------------------------
run(paste(sq(java_bin), java_mem, "-jar", sq(beagle_jar),
          paste0("gt=", sq(tmp)), paste0("out=", sq(out_prefix)),
          paste0("nthreads=", threads), paste0("seed=", CFG$seed)))
if (!file.exists(CFG$vcf_relaxed)) stop("Beagle finished but ", CFG$vcf_relaxed, " was not written - see ", out_prefix, ".log")
run(sprintf("%s index -f %s", sq(bcftools), sq(CFG$vcf_relaxed)))

# ---- 3. checks ------------------------------------------------------------------
n_out  <- trimws(system(sprintf("%s index -n %s", sq(bcftools), sq(CFG$vcf_relaxed)), intern = TRUE))
n_samp <- length(system(sprintf("%s query -l %s", sq(bcftools), sq(CFG$vcf_relaxed)), intern = TRUE))
bver   <- grep("^beagle", readLines(paste0(out_prefix, ".log"), n = 5), ignore.case = TRUE, value = TRUE)
msg("Relaxed, imputed SNPs:", n_out, "| samples:", n_samp, "| ", bver)
if (n_out != n_tmp) warning("SNP count changed during imputation (", n_tmp, " -> ", n_out, ")")
msg("Report this Beagle version in the Methods. Next: YAM_GENO_SET=relaxed Rscript 01_build_master_data.R")
