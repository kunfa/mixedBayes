#' Variable selection for a BayesQVGEL object
#'
#' Variable selection for a BayesQVGEL object
#'
#' @param obj BayesQVGEL object.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly..
#'
#' @details If sparse, the median probability model (MPM) (Barbieri and Berger, 2004) is used to identify predictors that are significantly associated
#' with the response variable. Otherwise, variable selection is based on 95\% credible interval.
#' Please check the references for more details about the variable selection.
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2022). Robust Bayesian variable selection for gene-environment interactions.
#' {\emph{Biometrics}, (in press)} \doi{10.1111/biom.13670}
#'
#' Barbieri, M.M. and Berger, J.O. (2004). Optimal predictive model selection. {\emph{Ann. Statist}, 32(3):870–897}
#'
#' @rdname selection
#' @return an object of class `selection' is returned, which is a list with component:
#' \item{inde}{a vector of indicators of selected effects.}
#'
#' @seealso \code{\link{BayesQVGEL}}
#'
#' @examples
#' data(data)
#' ## sparse
#' fit=BayesQVGEL(y,e,C,g,w,k,structure=c("group"))
#' selected=selection(fit,sparse=TRUE)
#' selected
#'
#' \donttest{
#' ## non-sparse
#' fit=BayesQVGEL(y,e,C,g,w,k,sparse=FALSE,structure=c("group"))
#' selected=selection(fit,sparse=FALSE)
#' selected
#' }
#'
#' @export
selection = function(obj,sparse){
  if(sparse){
    inde = selection_sparse(obj,burn.in=obj$burn.in)
  }
  else{
    inde = selection_nonsparse(obj,burn.in=obj$burn.in)
  }
  out = inde
  out
}
