#' fit a Bayesian quantile variable selection with mixed effects model in longitudinal study
#'
#' @keywords models
#' @param g the matrix of predictors (genetic factors) without intercept. Each row should be an observation vector.
#' @param y the matrix of response variable. The current version of rolong only supports continuous response.
#' @param e the matrix of a group of dummy environmental factors variables.
#' @param w the matrix of interactions between genetic factors and environmental factors.
#' @param iterations the number of MCMC iterations.
#' @param burn.in the number of iterations for burn-in.
#' @param slope logical flag. If TRUE, random intercept and slope model will be used.
#' @param robust logical flag. If TRUE, robust methods will be used.
#' @param quant specify different quantile when applying robust methods.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly.
#' @param structure structure for interaction effects, two choices are available. "group" for selection on group-level only. "individual" for selection on individual-level only.
#'
#' @details Consider the data model described in "\code{\link{data}}":
#' \deqn{Y_{ij} = X_{ij}^{T}\gamma_{0}+E_{ij}^{T}\gamma_{1}+\sum_{l=1}^{p}G_{ijl}\gamma_{2l}+\sum_{l=1}^{p}W_{ijl}^{T}\gamma_{3l}+Z_{ij}^{T}\alpha_{i}+\epsilon_{ij}.}
#' where \eqn{\gamma_{2l}} is the main effect of the \eqn{l}th genetic variant. The interaction effects is corresponding to the coefficient vector \eqn{\gamma_{3l}=(\gamma_{3l1}, \gamma_{3l2},\ldots,\gamma_{3lm})^\top}.
#'
#' When {structure="group"}, group-level selection will be conducted on \eqn{||\gamma_{3l}||_{2}}. If {structure="individual"}, individual-level selection will be conducted on each \eqn{\gamma_{3lq}}, (\eqn{q=1,\ldots,m}).
#'
#' When {slope=TRUE} (default), random intercept and slope model will be used as the mixed effects model.
#'
#' When {sparse=TRUE} (default), spike--and--slab priors are imposed on individual and/or group levels to identify important main and interaction effects. Otherwise, Laplacian shrinkage will be used.
#'
#' When {robust=TRUE} (default), the distribution of \eqn{\epsilon_{i}} is defined as a Laplace distribution with density.
#'
#' \eqn{
#' f(\epsilon_{ij}|\theta,\tau) = \theta(1-\theta)\exp\left\{-\tau\rho_{\theta}(\epsilon_{ij})\right\}
#' }, (\eqn{i=1,\dots,n,j=1,\dots,k }), which leads to a Bayesian formulation of quantile regression. If {robust=FALSE}, \eqn{\epsilon_{ij}} follows a normal distribution.
#' \cr
#'
#' Please check the references for more details about the prior distributions.
#'
#' @seealso \code{\link{data}}
#'
#' @examples
#' data(data)
#'
#' ## default method

#' fit=rolong(y,e,g,w,k, iterations=10000, burn.in=NULL, slope=TRUE, robust=TRUE, quant=0.5, sparse=TRUE, structure=c("group"))
#' fit$main
#' fit$inter
#'
#' \donttest{
#' ## alternative: robust individual selection
#' fit=rolong(y,e,g,w,k, iterations=10000, burn.in=NULL, slope=TRUE, robust=TRUE, quant=0.5, sparse=TRUE, structure=c("individual"))
#' fit$main
#' fit$inter
#'
#' ## alternative: non-robust group selection
#' fit=rolong(y,e,g,w,k, iterations=10000, burn.in=NULL, slope=TRUE, robust=FLASE, sparse=TRUE, structure=c("group"))
#' fit$main
#' fit$inter
#'
#' }
#'
#' @export

rolong <- function(y,e,g,w,k, iterations=10000, burn.in=NULL, slope=TRUE, robust=TRUE, quant=0.5, sparse=TRUE, structure=c("group","individual"))
{
  structure = match.arg(structure)
  this.call = match.call()

  if(iterations<1) stop("iterations must be a positive integer.")
  if(is.null(burn.in)){
    BI = floor(iterations/2)
  }else if(burn.in>=1){
    BI = as.integer(burn.in)
  }else{
    stop("burn.in must be a positive integer.")
  }
  if(iterations<=BI) stop("iterations must be larger than burn.in.")
  if(slope){
    z = cbind(rep(1,k),c(1:k))
  if(robust){
    out = LONRBGLSS(y,e,g,w,z,k,quant,iterations,sparse, structure)
  }else{
    out = LONBGLSS(y,e,g,w,z,k,iterations,sparse, structure)
  }
  }else{
    z = rep(1,k)
    if(robust){
      out = LONRBGLSS_1(y,e,g,w,z,k,quant,iterations,sparse, structure)
    }else{
      out = LONBGLSS_1(y,e,g,w,z,k,iterations,sparse, structure)
    }
  }

  coeff.main = apply(out$GS.beta[-(1:BI),,drop=FALSE], 2, stats::median);
  coeff.interaction = apply(out$GS.eta[-(1:BI),,drop=FALSE], 2, stats::median);


  coefficient = list(main=coeff.main, inter=coeff.interaction)

  fit = list(call = this.call, posterior = out, coefficient=coefficient, burn.in = BI, iterations=iterations)

  class(fit)=c("rolong", class(out))
  fit
}
