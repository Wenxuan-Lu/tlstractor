## usethis namespace: start
#' @useDynLib tlstractor, .registration = TRUE
#' @import data.table
#' @importFrom Rcpp evalCpp
## usethis namespace: end
#' @keywords internal
#' @noRd
"_PACKAGE"

# Package-private mutable state used internally (for example worker state on cluster nodes).
.tlstractor_env <- new.env(parent = emptyenv())

# Mark package as data.table-aware so DT[...] uses data.table semantics in package code.
.datatable.aware <- TRUE

utils::globalVariables(c(
	"AF", "ALT", "ALT_GDS", "BETA", "CHR", "CHR_ORIG", "GDS_ID", "ID", "POS", "REF", "REF_GDS", "SE",
	"i.ALT_GDS", "i.CHR_ORIG", "i.GDS_ID", "i.ID", "i.POS", "i.REF_GDS"
))

NULL
