#' This function changes the format of the longitudinal data from wide format to long format
#'
#' @param k the number of repeated measurement.
#' @param data either the response matrix or the predictor matrix.
#' @param type "r" for response matrix, "d" for design matrix.
#' @return the reformatted response vector or predictor matrix.
#' @usage reformat(k,data,type="r")
#' @export
reformat <- function(k, data, type = "r") {

  if (length(k) != 1 || !is.numeric(k) || is.na(k) || k <= 0 || k %% 1 != 0) {
    stop("k must be a positive integer.")
  }

  if (length(type) != 1 || !is.character(type) || is.na(type) ||
      !(type %in% c("r", "d"))) {
    stop("type must be exactly 'r' or 'd'.")
  }

  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  n <- nrow(data)

  if (type == "r") {

    if (ncol(data) != k) {
      stop("For type = 'r', data must be an n x k response matrix.")
    }

    out <- as.vector(t(data))

  } else if (type == "d") {

    if (ncol(data) == k) {
      stop("For type = 'd', data should not have ncol = k. Did you pass a response matrix?")
    }

    out <- data[rep(seq_len(n), each = k), , drop = FALSE]
  }

  return(out)
}
