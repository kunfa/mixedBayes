#' This function changes the format of the longitudinal data from wide format to long format
#'
#' @param k the number of repeated measurement.
#' @param data either the response matrix or the predictor matrix.
#' @param type "r" for response, "d" for design.
#' @return the reformatted response vector or predictor matrix.
#' @usage reformat(k,data,type=c("r","d"))
#' @export
reformat <- function(k,data,type=c("r","d")){
  type <- match.arg(type)
  if (!is.matrix(data)) data <- as.matrix(data)
  n <- nrow(data)
  if (type == "r") {

    if (ncol(data) != k)
      stop("For response type, ncol(data) must equal k.")

    return(as.vector(t(data)))
  }
  if (type == "d") {

    x <- data[rep(1:n, times = rep(k,n)), ]
    return(x)
  }

}
