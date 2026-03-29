#' Extract Local-Ancestry Tracts from Phased VCF and MSP (RFMix2/Gnomix)
#'
#' Parses a phased VCF together with an MSP local-ancestry file and writes
#' ancestry-specific outputs in one or more formats: `gds`, `txt`, `txt.gz`,
#' `vcf`, `vcf.gz`.
#'
#' @param vcf_path Character; input phased VCF path (`.vcf` or `.vcf.gz`).
#' @param msp_path Character; input MSP path.
#' @param num_ancs Integer; number of ancestral populations (must be > 1).
#' @param output_dir Character or `NULL`; output directory (must exist). If
#'   `NULL` (default), uses the directory of `vcf_path`.
#' @param output_formats Character vector; output format specification (for
#'   example `"gds"` or `c("gds", "txt.gz")`). Default is `"gds"`.
#'   Supported values: `gds`, `txt`, `txt.gz`, `vcf`, `vcf.gz`.
#'   If both compressed and uncompressed versions are requested for the same
#'   type, compressed output is written: `txt.gz` over `txt`, and `vcf.gz`
#'   over `vcf`.
#' @param chunk_size Integer; number of variants processed per chunk.
#'   Default is 1024; higher values can improve speed but increase memory usage.
#'
#' @return Invisibly returns `NULL`. Output files are written with filename
#'   prefix derived from `vcf_path` in `output_dir`.
#'
#' @details
#' Output files (prefix = basename of `vcf_path` without `.vcf` or `.vcf.gz`):
#' - `gds`: one file `<prefix>.gds` containing `sample.id`, variant metadata
#'   (`snp.chromosome`, `snp.position`, `snp.id`, `snp.ref`, `snp.alt`), and
#'   ancestry-specific nodes `dosage/anc0..ancK` and `hapcount/anc0..ancK`.
#' - `txt` / `txt.gz`: for each ancestry `k`, two files
#'   `<prefix>.anc{k}.dosage.txt(.gz)` and `<prefix>.anc{k}.hapcount.txt(.gz)`.
#'   Columns are `CHROM POS ID REF ALT` followed by one column per sample.
#' - `vcf` / `vcf.gz`: for each ancestry `k`, one file `<prefix>.anc{k}.vcf(.gz)`
#'   with `FORMAT=GT`; haplotypes not assigned to ancestry `k` are written as `.`.
#'
#' Assumptions for inputs:
#' - Biallelic variants.
#' - No missing genotype/local-ancestry fields.
#' - No duplicate variants.
#' - Autosomal chromosomes only (1-22, 01-22, `chr1`-`chr22`, `chr01`-`chr22`;
#'   case-insensitive).
#' - Alleles are uppercase `A/C/G/T`.
#' - `GT` is the first sample subfield (`GT:...`).
#' - `GT` is phased (`0|0`, `0|1`, `1|0`, `1|1`) with single-digit allele codes.
#' - One chromosome per VCF and per MSP file.
#' - Sample order is consistent between VCF and MSP files.
#'
#' @export
extract_tracts <- function(vcf_path, msp_path, num_ancs,
                           output_dir = NULL, output_formats = "gds", chunk_size = 1024L) {
  # Validate parameters
  validate_param_type(vcf_path, "character", "vcf_path", check_length_one = TRUE)

  if (!is.null(output_dir)) {
    validate_param_type(output_dir, "character", "output_dir", check_length_one = TRUE)
    if (!dir.exists(output_dir)) stop("Output directory does not exist: ", output_dir)
  }

  vinfo <- parse_vcf_filepath(vcf_path)
  out_dir <- if (is.null(output_dir)) vinfo$path else output_dir
  out_prefix <- file.path(out_dir, vinfo$prefix)

  validate_param_type(msp_path, "character", "msp_path", check_length_one = TRUE)
  if (!file.exists(msp_path)) stop("MSP file not found: ", msp_path)

  validate_param_type(num_ancs, "integer", "num_ancs", check_length_one = TRUE)
  num_ancs <- as.integer(num_ancs)
  if (num_ancs <= 1) stop("Error: num_ancs must be a positive integer greater than 1.", call. = FALSE)

  validate_param_type(output_formats, "character", "output_formats", check_length_one = FALSE)
  formats <- parse_tracts_output_formats(output_formats)

  validate_param_type(chunk_size, "integer", "chunk_size", check_length_one = TRUE)
  chunk_size <- as.integer(chunk_size)
  if (chunk_size <= 0) stop("Error: chunk_size must be a positive integer.", call. = FALSE)

  # Initialize GDS file if needed
  want_gds <- "gds" %in% formats
  gds_append_fn <- function(...) invisible(NULL)
  g <- NULL

  if (want_gds) {
    con <- open_in(vcf_path, vinfo$zipped)
    on.exit({
      try({
        if (inherits(con, "connection") && isOpen(con)) close(con)
      }, silent = TRUE)
    }, add = TRUE)
    sample_ids <- NULL
    repeat {
      ln <- readLines(con, n = 1L, warn = FALSE)
      if (!length(ln)) stop("EOF before VCF header.")
      if (startsWith(ln, "#CHROM") || (startsWith(ln, "#") && !startsWith(ln, "##"))) {
        if (endsWith(ln, "\t")) stop("Malformed VCF header: trailing tab.")
        hdr <- strsplit(ln, "\t", fixed = TRUE)[[1]]
        if (length(hdr) < 10L) stop("No samples found in VCF header.")
        sample_ids <- hdr[10:length(hdr)]
        break
      }
    }
    try({
      if (inherits(con, "connection") && isOpen(con)) close(con)
    }, silent = TRUE)

    gds_path <- sprintf("%s.gds", out_prefix)
    g <- init_gds(gds_path, sample_ids, num_ancs)
    on.exit({
      try({
        if (inherits(g, "gds.class")) {
          gdsfmt::closefn.gds(g)
          g <- NULL
        }
      }, silent = TRUE)
    }, add = TRUE)
    gds_append_fn <- make_gds_append_fn(g, num_ancs)
    rm(hdr, sample_ids)
  }

  # Call C++ function to process and write outputs
  extract_tracts_cpp(
    vcf_path = vcf_path,
    msp_path = msp_path,
    num_ancs = num_ancs,
    out_prefix = out_prefix,
    formats = formats,
    chunk_size = chunk_size,
    gds_append_fn = gds_append_fn
  )

  invisible(NULL)
}