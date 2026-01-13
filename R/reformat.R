#' This function changes the format of the longitudinal data from wide format to long format
#'
#' @param k the number of repeated measurement.
#' @param y the matrix of longitudinal response.
#' @param x the matrix of predictors.
#' @return a list containing the reformatted response vector and predictor matrix.
#' @usage reformat(k,y,x)
#' @export
reformat <- function(k,y,x){
  n=nrow(y)
  response=y
  y=rep(0,n*k)
  for (i in 1:n) {
    for (j in 1:k) {
      y[(i-1)*k+j]=response[i,j]
    }
  }

  data=cbind(y,x[rep(1:nrow(x), times = rep(k,n)), ])

  y=data[,1]
  x=data[,-1]
  return(list(y=y,x=x))
}
