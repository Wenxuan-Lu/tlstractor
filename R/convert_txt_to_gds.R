#' Convert Local-Ancestry Tracts Format from Text to GDS
#'
#' Converts ancestry-specific dosage/hapcount text files into a tracts GDS file.
#'
#' @param input_prefix Character; input file prefix used to locate ancestry files.
#'   Expected files per ancestry `k`: `<input_prefix>.anc{k}.dosage.txt(.gz)` and
#'   `<input_prefix>.anc{k}.hapcount.txt(.gz)`.
#' @param num_ancs Integer; number of ancestral populations (must be > 1).
#' @param output_prefix Character or `NULL`; output GDS prefix. If `NULL` (default),
#'   uses `input_prefix`.
#' @param chunk_size Integer; number of variants processed per chunk.
#'   Default is 1024; higher values can improve speed but increase memory usage.
#'
#' @return Invisibly returns `NULL`. Writes `<output_prefix>.gds`.
#'
#' @details
#' The generated GDS contains `sample.id`, variant metadata (`snp.chromosome`,
#' `snp.position`, `snp.id`, `snp.ref`, `snp.alt`), and ancestry-specific nodes
#' `dosage/anc0..ancK` and `hapcount/anc0..ancK`.
#'
#' @export
convert_txt_to_gds <- function(input_prefix, num_ancs, output_prefix=NULL, chunk_size = 1024L) {
  # Validate parameters
  validate_param_type(input_prefix, "character", "input_prefix", check_length_one = TRUE)

  validate_param_type(num_ancs, "integer", "num_ancs", check_length_one = TRUE)
  num_ancs <- as.integer(num_ancs)
  if (num_ancs <= 1) stop("Error: num_ancs must be a positive integer greater than 1.", call. = FALSE)

  if (!is.null(output_prefix)) {
    validate_param_type(output_prefix, "character", "output_prefix", check_length_one = TRUE)
    out_dir <- dirname(output_prefix)
    if (!dir.exists(out_dir)) stop("Output directory does not exist: ", out_dir)
  } else {
    output_prefix <- input_prefix
  }

  validate_param_type(chunk_size, "integer", "chunk_size", check_length_one = TRUE)
  chunk_size <- as.integer(chunk_size)
  if (chunk_size <= 0) stop("Error: chunk_size must be a positive integer.", call. = FALSE)

  # Build input file paths and validate
  files <- build_tracts_txt_input_paths(input_prefix, num_ancs)
  if (length(files$dosage_files) != length(files$hap_files)) {
    stop("Dosage/hapcount file count mismatch.", call. = FALSE)
  }
  message("Input dosage files: ", paste(files$dosage_files, collapse = ", "))
  message("Input hapcount files: ", paste(files$hap_files, collapse = ", "))
  sample_ids <- parse_tracts_txt_header(files$dosage_files[1L])
  if (length(sample_ids) <= 1L) {
    stop("Error: at least two samples need to be present.", call. = FALSE)
  }

  # Write to GDS file
  gds_path <- paste0(output_prefix, ".gds")
  g <- init_gds(gds_path, sample_ids, num_ancs)
  on.exit({
    try({
      if (inherits(g, "gds.class")) {
        gdsfmt::closefn.gds(g)
        g <- NULL
      }
    }, silent = TRUE)
  }, add = TRUE)

  nsamples <- length(sample_ids)
  rm(sample_ids)

  gds_append_fn <- make_gds_append_fn(g, num_ancs)
  convert_tracts_txt_to_gds_cpp(
    dosage_files = files$dosage_files,
    hap_files = files$hap_files,
    nsamples = nsamples,
    chunk_size = chunk_size,
    gds_append_fn = gds_append_fn
  )

  invisible(NULL)
}