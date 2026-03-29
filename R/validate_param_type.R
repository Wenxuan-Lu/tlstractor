#' @keywords internal
#' @noRd
validate_param_type <- function(value, expected_type, param_name, check_length_one = FALSE) {
  is_missing <- is.null(value) || anyNA(value) || (is.character(value) && any(value == ""))
  if (is_missing) {
    stop(sprintf("Error: Parameter '%s' is required and cannot be NULL, NA, or empty.", param_name), call. = FALSE)
  }
  if (check_length_one && length(value) != 1L) {
    stop(sprintf("Error: Parameter '%s' must have length 1, but has length %d.", param_name, length(value)), call. = FALSE)
  }

  if (expected_type %in% c("numeric", "integer")) {
    if (is.numeric(value) && any(!is.finite(value))) {
      stop(sprintf("Error: Parameter '%s' must be finite, but got non-finite values.", param_name), call. = FALSE)
    }
  }

  type_check <- switch(expected_type,
                       character = is.character(value),
                       integer = is.integer(value) || (is.numeric(value) && all(value == as.integer(value))),
                       numeric = is.numeric(value),
                       logical = is.logical(value),
                       FALSE
  )
  if (!type_check) {
    actual_type <- class(value)[1]
    stop(sprintf("Error: Parameter '%s' must be %s, but got %s.", param_name, expected_type, actual_type), call. = FALSE)
  }
}