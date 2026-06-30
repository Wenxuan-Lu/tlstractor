#' Munge GWAS Summary Statistics to Align with GDS File
#'
#' Aligns GWAS summary statistics with a GDS file using chromosome-position or variant ID matching.
#' Performs quality control filtering, allele flipping for discordant variants, and outputs
#' filtered results.
#'
#' @param gds_path Character; path to the input GDS file.
#' @param sumstats_path Character; path to the GWAS summary statistics file.
#' @param match_by Character; `"CHR-POS"` (default) or `"ID"`.
#'   `"CHR-POS"`: match by chromosome and position (autosomes 1-22 only).
#'   `"ID"`: match by variant ID.
#' @param chr_col Character; chromosome column name in the GWAS summary statistics file (default `"CHR"`).
#' @param pos_col Character; position column name in the GWAS summary statistics file (default `"POS"`).
#' @param id_col Character; variant ID column name in the GWAS summary statistics file (default `"ID"`).
#' @param ref_col Character; reference allele column name in the GWAS summary statistics file (default `"REF"`).
#' @param alt_col Character; alternative allele column name in the GWAS summary statistics file (default `"ALT"`).
#' @param beta_col Character; effect size column name in the GWAS summary statistics file (default `"BETA"`).
#' @param se_col Character; standard error column name in the GWAS summary statistics file (default `"SE"`).
#' @param n_col Character; total sample size column name in the GWAS summary statistics file (default `"N"`).
#'   Optional; omitted from output if absent in input GWAS summary statistics.
#' @param af_col Character; allele frequency column name in the GWAS summary statistics file (default `"AF"`).
#'   Optional; omitted from output if absent in input GWAS summary statistics.
#' @param n_case_col Character; case sample size column name in the GWAS summary statistics file (default `"N_case"`).
#'   Optional; omitted from output if absent in input GWAS summary statistics.
#' @param n_control_col Character; control sample size column name in the GWAS summary statistics file (default `"N_control"`).
#'   Optional; omitted from output if absent in input GWAS summary statistics.
#' @param af_case_col Character; case allele frequency column name in the GWAS summary statistics file (default `"AF_case"`).
#'   Optional; omitted from output if absent in input GWAS summary statistics.
#' @param af_control_col Character; control allele frequency column name in the GWAS summary statistics file (default `"AF_control"`).
#'   Optional; omitted from output if absent in input GWAS summary statistics.
#' @param remove_ambiguous Logical; if `TRUE` (default), exclude ambiguous SNPs.
#' @param output_path Character or `NULL`; output file path. If `NULL` (default),
#'   auto-generated as `<sumstats_prefix>_munged.txt` in the same directory as `sumstats_path`,
#'   where `<sumstats_prefix>` is the filename prefix of `sumstats_path` without extension.
#' @return Invisibly returns `NULL`. Writes filtered results to `output_path`
#'   as tab-delimited text. Output always includes `CHR`, `POS`, `ID`, `REF`, `ALT`, `BETA`, `SE`, and `GDS_ID`.
#'   Optional input columns `N`, `AF`, `N_case`, `N_control`, `AF_case`, and `AF_control` are preserved if present.
#'   `GDS_ID` is the 1-based variant index in the input GDS file to match variants
#'   in the summary statistics with variants in the GDS file.
#'
#' @details
#' **Required columns:**
#' The input summary statistics must contain columns for `REF`, `ALT`, `BETA`, and `SE`,
#' plus either (`CHR` and `POS`) when `match_by = "CHR-POS"` or `ID` when `match_by = "ID"`.
#' Optional columns are `N`, `AF`, `N_case`, `N_control`, `AF_case`, and `AF_control`.
#' Additional columns are allowed but ignored.
#' When `match_by = "CHR-POS"`, accepted chromosome formats are autosomes only:
#' `1-22`, `01-22`, `chr1-chr22`, or `chr01-chr22` (case-insensitive).
#' The `BETA` column should represent the effect size estimate for the alternative allele
#' (log odds ratio for logistic regression or linear coefficient for linear regression),
#' and the `SE` column should contain the corresponding standard error.
#'
#' **Quality control filters for summary statistics (applied in order):**
#'   1. REF/ALT: single A/C/G/T nucleotides, biallelic
#'   2. Ambiguous SNPs: optionally remove A/T, T/A, C/G, and G/C
#'   3. Autosomes: restrict to chromosomes 1-22 when `match_by = "CHR-POS"`
#'   4. POS: require `POS > 0` when `match_by = "CHR-POS"`
#'   5. ID: require non-missing, non-empty IDs when `match_by = "ID"`
#'   6. Duplicates: keep variants with unique `(CHR, POS)` or unique `ID` and remove all duplicates
#'   7. Effect/SE: require `BETA` not missing and `SE > 0`
#'   8. Sample sizes and AF ranges: require `N`, `N_case`, and `N_control` to be positive when present, and require `AF`, `AF_case`, and `AF_control` to lie in `(0, 1)` when present
#'
#' @export
munge_sumstats <- function(
    gds_path,
    sumstats_path,
    match_by = "CHR-POS",
    chr_col = "CHR",
    pos_col = "POS",
    id_col = "ID",
    ref_col = "REF",
    alt_col = "ALT",
    beta_col = "BETA",
    se_col = "SE",
    n_col = "N",
    af_col = "AF",
    n_case_col = "N_case",
    n_control_col = "N_control",
    af_case_col = "AF_case",
    af_control_col = "AF_control",
    remove_ambiguous = TRUE,
    output_path = NULL
) {

  # Validate inputs
  validate_param_type(gds_path, "character", "gds_path", check_length_one = TRUE)
  validate_param_type(sumstats_path, "character", "sumstats_path", check_length_one = TRUE)
  validate_param_type(match_by, "character", "match_by", check_length_one = TRUE)
  validate_param_type(chr_col, "character", "chr_col", check_length_one = TRUE)
  validate_param_type(pos_col, "character", "pos_col", check_length_one = TRUE)
  validate_param_type(id_col, "character", "id_col", check_length_one = TRUE)
  validate_param_type(ref_col, "character", "ref_col", check_length_one = TRUE)
  validate_param_type(alt_col, "character", "alt_col", check_length_one = TRUE)
  validate_param_type(beta_col, "character", "beta_col", check_length_one = TRUE)
  validate_param_type(se_col, "character", "se_col", check_length_one = TRUE)
  validate_param_type(n_col, "character", "n_col", check_length_one = TRUE)
  validate_param_type(af_col, "character", "af_col", check_length_one = TRUE)
  validate_param_type(n_case_col, "character", "n_case_col", check_length_one = TRUE)
  validate_param_type(n_control_col, "character", "n_control_col", check_length_one = TRUE)
  validate_param_type(af_case_col, "character", "af_case_col", check_length_one = TRUE)
  validate_param_type(af_control_col, "character", "af_control_col", check_length_one = TRUE)
  validate_param_type(remove_ambiguous, "logical", "remove_ambiguous", check_length_one = TRUE)
  if (is.null(output_path)) {
    output_path <- paste0(get_filepath_prefix(sumstats_path), "_munged.txt")
  } else {
    validate_param_type(output_path, "character", "output_path", check_length_one = TRUE)
  }
  if (!file.exists(gds_path)) {
    stop("Input GDS file not found: ", gds_path, call. = FALSE)
  }
  if (!file.exists(sumstats_path)) {
    stop("Summary statistics file not found: ", sumstats_path, call. = FALSE)
  }
  if (!(match_by %in% c("CHR-POS", "ID"))) {
    stop("match_by must be either 'CHR-POS' or 'ID'.", call. = FALSE)
  }

  # Validate summary statistics header
  sumstats_header <- data.table::fread(sumstats_path, nrows = 0L, data.table = FALSE)
  sumstats_names <- names(sumstats_header)
  if (is.null(sumstats_names) || length(sumstats_names) == 0L) {
    stop("Summary statistics file appears to have no header.", call. = FALSE)
  }

  if (match_by == "CHR-POS") {
    provided_cols <- c(chr_col, pos_col, ref_col, alt_col, beta_col, se_col,
                       n_col, af_col, n_case_col, n_control_col, af_case_col, af_control_col)
    required_cols <- c(chr_col, pos_col, ref_col, alt_col, beta_col, se_col)
    new_names <- c("CHR", "POS", "REF", "ALT", "BETA", "SE")
  } else {
    provided_cols <- c(id_col, ref_col, alt_col, beta_col, se_col,
                       n_col, af_col, n_case_col, n_control_col, af_case_col, af_control_col)
    required_cols <- c(id_col, ref_col, alt_col, beta_col, se_col)
    new_names <- c("ID", "REF", "ALT", "BETA", "SE")
  }

  invalid_cols <- provided_cols[make.names(provided_cols) != provided_cols]
  if (length(invalid_cols) > 0L) {
    stop(
      "Invalid column name(s): ", paste(invalid_cols, collapse = ", "),
      ". Use valid R column names (letters, numbers, dots, or underscores; no leading digit).",
      call. = FALSE
    )
  }

  missing_cols <- setdiff(required_cols, sumstats_names)
  if (length(missing_cols) > 0L) {
    stop("Missing required summary statistics column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  optional_cols <- c(n_col, af_col, n_case_col, n_control_col, af_case_col, af_control_col)
  optional_out_names <- c("N", "AF", "N_case", "N_control", "AF_case", "AF_control")
  optional_present <- optional_cols %in% sumstats_names
  selected_cols <- c(required_cols, optional_cols[optional_present])
  new_names <- c(new_names, optional_out_names[optional_present])

  # Set column classes for efficient reading
  numeric_optional_cols <- optional_cols[optional_present]
  if (match_by == "CHR-POS") {
    colClasses <- list(
      character = c(chr_col, ref_col, alt_col),
      integer   = pos_col,
      numeric   = c(beta_col, se_col, numeric_optional_cols)
    )
  } else {
    colClasses <- list(
      character = c(id_col, ref_col, alt_col),
      numeric   = c(beta_col, se_col, numeric_optional_cols)
    )
  }

  # Read summary statistics
  sumstats <- data.table::fread(
    sumstats_path,
    select = selected_cols,
    colClasses = colClasses,
    data.table = TRUE,
    showProgress = FALSE
  )
  data.table::setnames(sumstats, old = selected_cols, new = new_names, skip_absent = FALSE)
  data.table::setcolorder(sumstats, new_names)
  if (nrow(sumstats) == 0L) {
    stop("No rows read from summary statistics file.", call. = FALSE)
  }

  # Read GDS file
  in_gds <- gdsfmt::openfn.gds(gds_path, readonly = TRUE)
  on.exit({
    try({
      if (inherits(in_gds, "gds.class")) {
        gdsfmt::closefn.gds(in_gds)
        in_gds <- NULL
      }
    }, silent = TRUE)
  }, add = TRUE)

  chrom_node <- gdsfmt::index.gdsn(in_gds, "snp.chromosome")
  pos_node <- gdsfmt::index.gdsn(in_gds, "snp.position")
  id_node <- gdsfmt::index.gdsn(in_gds, "snp.id")
  ref_node <- gdsfmt::index.gdsn(in_gds, "snp.ref")
  alt_node <- gdsfmt::index.gdsn(in_gds, "snp.alt")

  gds_chr_orig <- gdsfmt::read.gdsn(chrom_node)

  gds_dt <- data.table::data.table(
    CHR_ORIG = gds_chr_orig,
    POS = gdsfmt::read.gdsn(pos_node),
    ID = gdsfmt::read.gdsn(id_node),
    REF = toupper(gdsfmt::read.gdsn(ref_node)),
    ALT = toupper(gdsfmt::read.gdsn(alt_node)),
    GDS_ID = seq_along(gds_chr_orig)
  )

  if (match_by == "CHR-POS") {
    gds_chr_clean <- sub("^[Cc][Hh][Rr]", "", gds_chr_orig)
    gds_chr_int <- suppressWarnings(as.integer(gds_chr_clean))
    data.table::set(gds_dt, j = "CHR", value = gds_chr_int)
    rm(gds_chr_clean, gds_chr_int)
  }
  rm(gds_chr_orig)

  if (nrow(gds_dt) == 0L) {
    stop("No variants found in the GDS file.", call. = FALSE)
  }

  # Perform QC and matching
  if (match_by == "CHR-POS") {
    out <- qc_and_match_by_chr_pos(
      sumstats,
      gds_dt,
      optional_out_names[optional_present],
      remove_ambiguous
    )
  } else {
    out <- qc_and_match_by_id(
      sumstats,
      gds_dt,
      optional_out_names[optional_present],
      remove_ambiguous
    )
  }
  rm(gds_dt, sumstats)

  if (nrow(out) == 0L) {
    stop("No overlapping variants between GDS and summary statistics after filtering.", call. = FALSE)
  }

  # Finalize output
  data.table::setnames(out, "CHR_ORIG", "CHR")
  data.table::setcolorder(out, c("CHR", "POS", "ID", "REF", "ALT", "BETA", "SE", optional_out_names[optional_present], "GDS_ID"))
  data.table::setorder(out, GDS_ID)
  data.table::fwrite(out, file = output_path, sep = "\t", quote = FALSE, col.names = TRUE)
  message("Variants retained after QC and matching: ", nrow(out))
  message("Munged summary statistics written to: ", output_path)

  invisible(NULL)
}


#' Get prefix of a fread-readable file path
#'
#' Removes compression (e.g., .gz, .bz2, .xz, .zip) and one file extension.
#'
#' @keywords internal
#' @noRd
get_filepath_prefix <- function(filepath) {
  base <- basename(filepath)
  dir  <- dirname(filepath)
  base <- sub("\\.(gz|bz2|xz|zip)$", "", base, ignore.case = TRUE)
  base <- sub("(?<!^)\\.[^.]+$", "", base, perl = TRUE)
  file.path(dir, base)
}


#' @keywords internal
#' @noRd
qc_and_match_by_chr_pos <- function(sumstats, gds_dt,
                                    optional_cols,
                                    remove_ambiguous) {
  # Validate REF/ALT
  data.table::set(sumstats, j = "REF", value = toupper(sumstats[["REF"]]))
  data.table::set(sumstats, j = "ALT", value = toupper(sumstats[["ALT"]]))
  valid_alleles <- c("A", "C", "G", "T")
  keep <- (nchar(sumstats[["REF"]]) == 1L) &
    (nchar(sumstats[["ALT"]]) == 1L) &
    (data.table::`%chin%`(sumstats[["REF"]], valid_alleles)) &
    (data.table::`%chin%`(sumstats[["ALT"]], valid_alleles)) &
    (sumstats[["REF"]] != sumstats[["ALT"]])
  invalid_count <- sum(!keep)
  message(invalid_count, " rows excluded: invalid alleles (not biallelic SNPs).")
  if (invalid_count > 0) sumstats <- sumstats[keep]

  # Remove ambiguous SNPs
  if (remove_ambiguous) {
    ref_alt <- paste0(sumstats[["REF"]], sumstats[["ALT"]])
    keep <- !(data.table::`%chin%`(ref_alt, c("AT", "TA", "CG", "GC")))
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: ambiguous SNPs (A/T, T/A, C/G, G/C).")
    if (invalid_count > 0) sumstats <- sumstats[keep]
    rm(ref_alt)
  }

  # Process chromosomes: remove chr/CHR prefix, convert to integer, keep autosomes 1-22
  chr_vec <- sub("^[Cc][Hh][Rr]", "", sumstats[["CHR"]])
  chr_int <- suppressWarnings(as.integer(chr_vec))
  keep <- !is.na(chr_int) & chr_int >= 1L & chr_int <= 22L
  invalid_count <- sum(!keep)
  message(invalid_count, " rows excluded: non-autosomal or invalid chromosomes.")
  if (invalid_count > 0) {
    sumstats <- sumstats[keep]
    chr_int <- chr_int[keep]
  }
  data.table::set(sumstats, j = "CHR", value = chr_int)
  rm(chr_vec, chr_int)

  # Validate POS
  keep <- !is.na(sumstats[["POS"]]) & sumstats[["POS"]] > 0L
  invalid_count <- sum(!keep)
  message(invalid_count, " rows excluded: invalid POS values.")
  if (invalid_count > 0) sumstats <- sumstats[keep]

  # Remove duplicates
  n_before <- nrow(sumstats)
  sumstats <- sumstats[
    !duplicated(sumstats, by = c("CHR", "POS")) &
    !duplicated(sumstats, by = c("CHR", "POS"), fromLast = TRUE)
  ]
  message(n_before - nrow(sumstats), " rows excluded: duplicate (CHR, POS) pairs.")

  # Validate BETA, SE
  keep <- !is.na(sumstats[["BETA"]]) & !is.na(sumstats[["SE"]]) & sumstats[["SE"]] > 0
  invalid_count <- sum(!keep)
  message(invalid_count, " rows excluded: invalid BETA or SE values.")
  if (invalid_count > 0) sumstats <- sumstats[keep]

  # Validate N, AF, N_case, N_control, AF_case, AF_control if present
  if ("N" %in% optional_cols) {
    keep <- !is.na(sumstats[["N"]]) & sumstats[["N"]] > 0
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid N values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("AF" %in% optional_cols) {
    keep <- !is.na(sumstats[["AF"]]) & sumstats[["AF"]] > 0 & sumstats[["AF"]] < 1
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid AF values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("N_case" %in% optional_cols) {
    keep <- !is.na(sumstats[["N_case"]]) & sumstats[["N_case"]] > 0
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid N_case values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("N_control" %in% optional_cols) {
    keep <- !is.na(sumstats[["N_control"]]) & sumstats[["N_control"]] > 0
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid N_control values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("AF_case" %in% optional_cols) {
    keep <- !is.na(sumstats[["AF_case"]]) & sumstats[["AF_case"]] > 0 & sumstats[["AF_case"]] < 1
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid AF_case values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("AF_control" %in% optional_cols) {
    keep <- !is.na(sumstats[["AF_control"]]) & sumstats[["AF_control"]] > 0 & sumstats[["AF_control"]] < 1
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid AF_control values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  rm(keep)

  if (nrow(sumstats) == 0L) {
    stop("No rows remain in summary statistics after filtering.", call. = FALSE)
  }

  # Index sumstats for faster joining
  data.table::setindexv(sumstats, c("CHR", "POS", "REF", "ALT"))

  # Forward match
  fwd_match <- sumstats[gds_dt, on = c("CHR", "POS", "REF", "ALT"), nomatch = 0,
    c(
      list(
        CHR_ORIG = i.CHR_ORIG,
        POS      = i.POS,
        ID       = i.ID,
        REF      = i.REF,
        ALT      = i.ALT
      ),
      mget(c("BETA", "SE", optional_cols)),
      list(
        GDS_ID = i.GDS_ID
      )
    )
  ]
  
  # Track matched GDS IDs
  matched <- logical(nrow(gds_dt))
  if (nrow(fwd_match) > 0L) {
    matched[fwd_match[["GDS_ID"]]] <- TRUE
  }

  # Flip match for unmatched GDS variants
  if (any(!matched)) {
    flip_match <- sumstats[gds_dt[!matched], on = c("CHR", "POS", "REF" = "ALT", "ALT" = "REF"), nomatch = 0,
      c(
        list(
          CHR_ORIG = i.CHR_ORIG,
          POS      = i.POS,
          ID       = i.ID,
          REF      = i.REF,
          ALT      = i.ALT
        ),
        mget(c("BETA", "SE", optional_cols)),
        list(
          GDS_ID = i.GDS_ID
        )
      )
    ]
    data.table::set(flip_match, j = "BETA", value = -flip_match[["BETA"]])
    if ("AF" %in% optional_cols) {
      data.table::set(flip_match, j = "AF", value = 1 - flip_match[["AF"]])
    }
    if ("AF_case" %in% optional_cols) {
      data.table::set(flip_match, j = "AF_case", value = 1 - flip_match[["AF_case"]])
    }
    if ("AF_control" %in% optional_cols) {
      data.table::set(flip_match, j = "AF_control", value = 1 - flip_match[["AF_control"]])
    }
    out <- data.table::rbindlist(list(fwd_match, flip_match), use.names = TRUE)
  } else {
    out <- fwd_match
  }

  return(out)
}


#' @keywords internal
#' @noRd
qc_and_match_by_id <- function(sumstats, gds_dt,
                               optional_cols,
                               remove_ambiguous) {
  # Validate REF/ALT
  data.table::set(sumstats, j = "REF", value = toupper(sumstats[["REF"]]))
  data.table::set(sumstats, j = "ALT", value = toupper(sumstats[["ALT"]]))
  valid_alleles <- c("A", "C", "G", "T")
  keep <- (nchar(sumstats[["REF"]]) == 1L) &
    (nchar(sumstats[["ALT"]]) == 1L) &
    (data.table::`%chin%`(sumstats[["REF"]], valid_alleles)) &
    (data.table::`%chin%`(sumstats[["ALT"]], valid_alleles)) &
    (sumstats[["REF"]] != sumstats[["ALT"]])
  invalid_count <- sum(!keep)
  message(invalid_count, " rows excluded: invalid alleles (not biallelic SNPs).")
  if (invalid_count > 0) sumstats <- sumstats[keep]

  # Remove ambiguous SNPs
  if (remove_ambiguous) {
    ref_alt <- paste0(sumstats[["REF"]], sumstats[["ALT"]])
    keep <- !(data.table::`%chin%`(ref_alt, c("AT", "TA", "CG", "GC")))
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: ambiguous SNPs (A/T, T/A, C/G, G/C).")
    if (invalid_count > 0) sumstats <- sumstats[keep]
    rm(ref_alt)
  }

  # Validate ID column
  keep <- !is.na(sumstats[["ID"]]) & sumstats[["ID"]] != ""
  invalid_count <- sum(!keep)
  message(invalid_count, " rows excluded: missing or empty ID values.")
  if (invalid_count > 0) sumstats <- sumstats[keep]

  # Remove duplicates
  n_before <- nrow(sumstats)
  sumstats <- sumstats[
    !duplicated(sumstats[["ID"]]) &
    !duplicated(sumstats[["ID"]], fromLast = TRUE)
  ]
  message(n_before - nrow(sumstats), " rows excluded: duplicate IDs.")

  # Validate BETA, SE
  keep <- !is.na(sumstats[["BETA"]]) & !is.na(sumstats[["SE"]]) & sumstats[["SE"]] > 0
  invalid_count <- sum(!keep)
  message(invalid_count, " rows excluded: invalid BETA or SE values.")
  if (invalid_count > 0) sumstats <- sumstats[keep]

  # Validate N, AF, N_case, N_control, AF_case, AF_control if present
  if ("N" %in% optional_cols) {
    keep <- !is.na(sumstats[["N"]]) & sumstats[["N"]] > 0
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid N values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("AF" %in% optional_cols) {
    keep <- !is.na(sumstats[["AF"]]) & sumstats[["AF"]] > 0 & sumstats[["AF"]] < 1
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid AF values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("N_case" %in% optional_cols) {
    keep <- !is.na(sumstats[["N_case"]]) & sumstats[["N_case"]] > 0
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid N_case values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("N_control" %in% optional_cols) {
    keep <- !is.na(sumstats[["N_control"]]) & sumstats[["N_control"]] > 0
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid N_control values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("AF_case" %in% optional_cols) {
    keep <- !is.na(sumstats[["AF_case"]]) & sumstats[["AF_case"]] > 0 & sumstats[["AF_case"]] < 1
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid AF_case values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  if ("AF_control" %in% optional_cols) {
    keep <- !is.na(sumstats[["AF_control"]]) & sumstats[["AF_control"]] > 0 & sumstats[["AF_control"]] < 1
    invalid_count <- sum(!keep)
    message(invalid_count, " rows excluded: invalid AF_control values.")
    if (invalid_count > 0) sumstats <- sumstats[keep]
  }
  rm(keep)

  if (nrow(sumstats) == 0L) {
    stop("No rows remain in summary statistics after filtering.", call. = FALSE)
  }

  # Index sumstats for faster joining
  data.table::setindexv(sumstats, c("ID", "REF", "ALT"))

  # Forward match: CHR and POS come from GDS
  fwd_match <- sumstats[gds_dt, on = c("ID", "REF", "ALT"), nomatch = 0,
    c(
      list(
        CHR_ORIG = i.CHR_ORIG,
        POS      = i.POS,
        ID       = i.ID,
        REF      = i.REF,
        ALT      = i.ALT
      ),
      mget(c("BETA", "SE", optional_cols)),
      list(
        GDS_ID = i.GDS_ID
      )
    )
  ]

  # Track matched GDS IDs
  matched <- logical(nrow(gds_dt))
  if (nrow(fwd_match) > 0L) {
    matched[fwd_match[["GDS_ID"]]] <- TRUE
  }

  # Flip match for unmatched GDS variants
  if (any(!matched)) {
    flip_match <- sumstats[gds_dt[!matched], on = c("ID", "REF" = "ALT", "ALT" = "REF"), nomatch = 0,
      c(
        list(
          CHR_ORIG = i.CHR_ORIG,
          POS      = i.POS,
          ID       = i.ID,
          REF      = i.REF,
          ALT      = i.ALT
        ),
        mget(c("BETA", "SE", optional_cols)),
        list(
          GDS_ID = i.GDS_ID
        )
      )
    ]
    data.table::set(flip_match, j = "BETA", value = -flip_match[["BETA"]])
    if ("AF" %in% optional_cols) {
      data.table::set(flip_match, j = "AF", value = 1 - flip_match[["AF"]])
    }
    if ("AF_case" %in% optional_cols) {
      data.table::set(flip_match, j = "AF_case", value = 1 - flip_match[["AF_case"]])
    }
    if ("AF_control" %in% optional_cols) {
      data.table::set(flip_match, j = "AF_control", value = 1 - flip_match[["AF_control"]])
    }
    out <- data.table::rbindlist(list(fwd_match, flip_match), use.names = TRUE)
  } else {
    out <- fwd_match
  }

  return(out)
}
