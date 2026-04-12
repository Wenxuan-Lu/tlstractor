#' Convert Local-Ancestry Tracts Format from GDS to Text
#'
#' Converts ancestry-specific dosage and hapcount matrices stored in a tracts
#' GDS file into tab-delimited text files (plain or gzipped).
#'
#' @param gds_path Character; input tracts GDS path.
#' @param output_prefix Character or `NULL`; output file prefix. If `NULL` (default),
#'   uses the input filename prefix in the same directory as `gds_path`.
#' @param output_format Character; output format, one of `"txt"` or `"txt.gz"` (default).
#' @param chunk_size Integer; number of variants processed per chunk.
#'   Default is 1024; higher values can improve speed but increase memory usage.
#'
#' @return Invisibly returns `NULL`. For each ancestry `k`, writes two files:
#'   `<output_prefix>.anc{k}.dosage.txt(.gz)` and
#'   `<output_prefix>.anc{k}.hapcount.txt(.gz)`.
#'
#' @details
#' Output tables contain columns `CHROM`, `POS`, `ID`, `REF`, `ALT`, followed by
#' one column per sample from `sample.id` in the input GDS.
#'
#' @export
convert_gds_to_txt <- function(gds_path, output_prefix=NULL, output_format="txt.gz", chunk_size = 1024L) {
  # Validate parameters
  validate_param_type(gds_path, "character", "gds_path", check_length_one = TRUE)
  if (!file.exists(gds_path)) stop("Input GDS file not found: ", gds_path)

  if (!is.null(output_prefix)) {
    validate_param_type(output_prefix, "character", "output_prefix", check_length_one = TRUE)
    out_dir <- dirname(output_prefix)
    if (!dir.exists(out_dir)) stop("Output directory does not exist: ", out_dir)
  } else {
    output_prefix <- file.path(dirname(gds_path), tools::file_path_sans_ext(basename(gds_path)))
  }

  validate_param_type(output_format, "character", "output_format", check_length_one = TRUE)
  fmt <- match.arg(output_format, c("txt", "txt.gz"))
  gz <- fmt == "txt.gz"

  validate_param_type(chunk_size, "integer", "chunk_size", check_length_one = TRUE)
  chunk_size <- as.integer(chunk_size)
  if (chunk_size <= 0) stop("Error: chunk_size must be a positive integer.", call. = FALSE)

  # Open GDS file and read metadata
  in_gds <- gdsfmt::openfn.gds(gds_path, readonly = TRUE)
  on.exit({
    try({
      if (inherits(in_gds, "gds.class")) {
        gdsfmt::closefn.gds(in_gds)
        in_gds <- NULL
      }
    }, silent = TRUE)
  }, add = TRUE)

  sample_ids <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(in_gds, "sample.id"))
  nsamples <- length(sample_ids)
  if (is.null(nsamples) || is.na(nsamples) || nsamples <= 1L) {
    stop("Error: at least two samples need to be present.", call. = FALSE)
  }

  chrom_node <- gdsfmt::index.gdsn(in_gds, "snp.chromosome")
  pos_node <- gdsfmt::index.gdsn(in_gds, "snp.position")
  id_node <- gdsfmt::index.gdsn(in_gds, "snp.id")
  ref_node <- gdsfmt::index.gdsn(in_gds, "snp.ref")
  alt_node <- gdsfmt::index.gdsn(in_gds, "snp.alt")

  num_ancs <- length(gdsfmt::ls.gdsn(gdsfmt::index.gdsn(in_gds, "dosage")))
  if (is.null(num_ancs) || is.na(num_ancs) || num_ancs <= 1L) {
    stop("Error: at least two ancestries need to be present under 'dosage' in GDS file.", call. = FALSE)
  }
  nvars <- gdsfmt::objdesp.gdsn(chrom_node)$dim
  if (is.null(nvars)) nvars <- gdsfmt::objdesp.gdsn(chrom_node)$size
  nvars <- as.integer(nvars)

  # Initialize output writers and write headers
  dosage_paths <- character(num_ancs)
  hapcount_paths <- character(num_ancs)
  dosage_writers <- vector("list", num_ancs)
  hapcount_writers <- vector("list", num_ancs)
  on.exit({
    for (k in seq_len(num_ancs)) {
      if (!is.null(dosage_writers[[k]])) gds_to_txt_close_writer_cpp(dosage_writers[[k]])
    }
    for (k in seq_len(num_ancs)) {
      if (!is.null(hapcount_writers[[k]])) gds_to_txt_close_writer_cpp(hapcount_writers[[k]])
    }
  }, add = TRUE)
  for (k in seq_len(num_ancs)) {
    dosage_paths[k] <- sprintf("%s.anc%d.dosage.txt%s", output_prefix, k - 1, if (gz) ".gz" else "")
    hapcount_paths[k] <- sprintf("%s.anc%d.hapcount.txt%s", output_prefix, k - 1, if (gz) ".gz" else "")
    dosage_writers[[k]] <- gds_to_txt_open_writer_cpp(dosage_paths[k], gz)
    hapcount_writers[[k]] <- gds_to_txt_open_writer_cpp(hapcount_paths[k], gz)
    gds_to_txt_write_header_cpp(dosage_writers[[k]], sample_ids)
    gds_to_txt_write_header_cpp(hapcount_writers[[k]], sample_ids)
  }
  rm(sample_ids)

  # Get GDS nodes for dosage and hapcount matrices
  dosage_nodes <- lapply(0:(num_ancs - 1L), function(k) gdsfmt::index.gdsn(in_gds, paste0("dosage/anc", k)))
  hapcount_nodes <- lapply(0:(num_ancs - 1L), function(k) gdsfmt::index.gdsn(in_gds, paste0("hapcount/anc", k)))

  # Process variants in chunks and write to text files
  for (start in seq(1, nvars, by = chunk_size)) {
    end <- min(start + chunk_size - 1, nvars)
    idx_len <- end - start + 1L
    chroms <- gdsfmt::read.gdsn(chrom_node, start = start, count = idx_len)
    pos <- gdsfmt::read.gdsn(pos_node, start = start, count = idx_len)
    ids <- gdsfmt::read.gdsn(id_node, start = start, count = idx_len)
    refs <- gdsfmt::read.gdsn(ref_node, start = start, count = idx_len)
    alts <- gdsfmt::read.gdsn(alt_node, start = start, count = idx_len)

    for (k in seq_len(num_ancs)) {
      mat <- gdsfmt::read.gdsn(dosage_nodes[[k]], start = c(1, start), count = c(nsamples, idx_len))
      if (!is.matrix(mat)) {
        dim(mat) <- c(nsamples, idx_len)
      }
      gds_to_txt_append_chunk_cpp(
        writer_ptr = dosage_writers[[k]],
        chroms = chroms,
        pos = pos,
        ids = ids,
        refs = refs,
        alts = alts,
        mat = mat
      )

      mat <- gdsfmt::read.gdsn(hapcount_nodes[[k]], start = c(1, start), count = c(nsamples, idx_len))
      if (!is.matrix(mat)) {
        dim(mat) <- c(nsamples, idx_len)
      }
      gds_to_txt_append_chunk_cpp(
        writer_ptr = hapcount_writers[[k]],
        chroms = chroms,
        pos = pos,
        ids = ids,
        refs = refs,
        alts = alts,
        mat = mat
      )
    }
  }

  invisible(NULL)
}