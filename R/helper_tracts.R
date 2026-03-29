#' @keywords internal
#' @noRd
open_in <- function(path, zipped) if (zipped) gzfile(path, "rt") else file(path, "rt")


#' @keywords internal
#' @noRd
init_gds <- function(gds_path, sample_ids, num_ancs) {
  g <- gdsfmt::createfn.gds(gds_path)
  gdsfmt::add.gdsn(g, "sample.id", sample_ids, compress = "ZIP_RA")

  gdsfmt::add.gdsn(g, "snp.chromosome", character(), storage = "string", compress = "ZIP_RA")
  gdsfmt::add.gdsn(g, "snp.position", integer(), storage = "int32", compress = "ZIP_RA")
  gdsfmt::add.gdsn(g, "snp.id", character(), storage = "string", compress = "ZIP_RA")
  gdsfmt::add.gdsn(g, "snp.ref", character(), storage = "string", compress = "ZIP_RA")
  gdsfmt::add.gdsn(g, "snp.alt", character(), storage = "string", compress = "ZIP_RA")

  dos_folder <- gdsfmt::addfolder.gdsn(g, "dosage")
  hap_folder <- gdsfmt::addfolder.gdsn(g, "hapcount")

  ns <- length(sample_ids)
  for (k in 0:(num_ancs - 1L)) {
    gdsfmt::add.gdsn(dos_folder, paste0("anc", k),
             val = array(as.raw(integer()), dim = c(ns, 0L)),
             storage = "bit2", compress = "LZ4_RA.hc:256K")
    gdsfmt::add.gdsn(hap_folder, paste0("anc", k),
             val = array(as.raw(integer()), dim = c(ns, 0L)),
             storage = "bit2", compress = "LZ4_RA.hc:256K")
  }
  g
}


#' @keywords internal
#' @noRd
make_gds_append_fn <- function(g, num_ancs) {
  force(g)
  chrom_node <- gdsfmt::index.gdsn(g, "snp.chromosome")
  pos_node <- gdsfmt::index.gdsn(g, "snp.position")
  id_node <- gdsfmt::index.gdsn(g, "snp.id")
  ref_node <- gdsfmt::index.gdsn(g, "snp.ref")
  alt_node <- gdsfmt::index.gdsn(g, "snp.alt")
  dos_nodes <- lapply(0:(num_ancs - 1L), function(k) gdsfmt::index.gdsn(g, paste0("dosage/anc", k)))
  hap_nodes <- lapply(0:(num_ancs - 1L), function(k) gdsfmt::index.gdsn(g, paste0("hapcount/anc", k)))
  function(chrom, pos, id, ref, alt, dosage_list, hap_list) {
    gdsfmt::append.gdsn(chrom_node, chrom)
    gdsfmt::append.gdsn(pos_node, pos)
    gdsfmt::append.gdsn(id_node, id)
    gdsfmt::append.gdsn(ref_node, ref)
    gdsfmt::append.gdsn(alt_node, alt)

    for (k in 0:(num_ancs - 1L)) {
      gdsfmt::append.gdsn(dos_nodes[[k + 1L]], dosage_list[[k + 1L]])
      gdsfmt::append.gdsn(hap_nodes[[k + 1L]], hap_list[[k + 1L]])
    }
    invisible(NULL)
  }
}


#' @keywords internal
#' @noRd
parse_vcf_filepath <- function(vcf_path) {
  if (!file.exists(vcf_path)) stop("File not found: ", vcf_path)
  if (dir.exists(vcf_path)) stop("Expected a file, but found a directory: ", vcf_path)
  filename <- basename(vcf_path)
  path <- dirname(vcf_path)

  if (endsWith(filename, ".vcf.gz")) {
    list(path = path, prefix = sub("\\.vcf\\.gz$", "", filename), zipped = TRUE)
  } else if (endsWith(filename, ".vcf")) {
    list(path = path, prefix = sub("\\.vcf$", "", filename), zipped = FALSE)
  } else {
    stop("Unexpected file extension: ", filename, ". VCF file must end in .vcf or .vcf.gz")
  }
}


#' @keywords internal
#' @noRd
parse_tracts_output_formats <- function(fmt) {
  supported_formats <- c("gds", "txt", "txt.gz", "vcf", "vcf.gz")

  if (is.null(fmt) || !is.character(fmt) || anyNA(fmt) || any(trimws(fmt) == "")) {
    stop(
      "output_formats must be a character vector with supported formats: ",
      paste(supported_formats, collapse = ", "),
      call. = FALSE
    )
  }

  fmt <- trimws(fmt)

  invalid <- setdiff(fmt, supported_formats)
  if (length(invalid) > 0) {
    stop(
      "Unsupported output format(s): ",
      paste(invalid, collapse = ", "),
      "\nSupported formats are: ",
      paste(supported_formats, collapse = ", "),
      call. = FALSE
    )
  }

  unique(fmt)
}


#' @keywords internal
#' @noRd
build_tracts_txt_input_paths <- function(prefix, num_ancs) {
  dosage_files <- character(num_ancs)
  hap_files <- character(num_ancs)
  for (k in 0:(num_ancs - 1L)) {
    base_d <- sprintf("%s.anc%d.dosage.txt", prefix, k)
    base_h <- sprintf("%s.anc%d.hapcount.txt", prefix, k)
    d <- if (file.exists(paste0(base_d, ".gz"))) paste0(base_d, ".gz") else base_d
    h <- if (file.exists(paste0(base_h, ".gz"))) paste0(base_h, ".gz") else base_h
    if (!file.exists(d)) stop("Missing dosage file for ancestry ", k, call. = FALSE)
    if (!file.exists(h)) stop("Missing hapcount file for ancestry ", k, call. = FALSE)
    dosage_files[k + 1L] <- d
    hap_files[k + 1L] <- h
  }
  list(dosage_files = dosage_files, hap_files = hap_files)
}


#' @keywords internal
#' @noRd
parse_tracts_txt_header <- function(path) {
  con <- open_in(path, endsWith(path, ".gz"))
  on.exit({
    try({
      if (inherits(con, "connection") && isOpen(con)) close(con)
    }, silent = TRUE)
  }, add = TRUE)
  ln <- readLines(con, n = 1L, warn = FALSE)
  if (!length(ln)) stop("Empty txt file: ", path, call. = FALSE)
  if (endsWith(ln, "\t")) stop("Malformed txt header: trailing tab.", call. = FALSE)
  hdr <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  if (length(hdr) < 6L) stop("No samples found in txt header.", call. = FALSE)
  hdr[6:length(hdr)]
}