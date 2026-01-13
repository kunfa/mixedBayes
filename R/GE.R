#' Construct Gene-Environment (G×E) Interaction Matrix
#' @param g the long-format matrix of genetic predictors.
#' @param e the long-format design matrix for environment/treatment effects.
#' @return the G×E interaction terms.
#' @usage GE(g,e)
#' @export
#'
GE = function(g,e){
  m = ncol(g)
  w = c()

  for (i in 1:m)
  {

    w = cbind(w,g[,i]*e)

  }
  return(w)
}
