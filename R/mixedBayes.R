#' fit a Bayesian longitudinal regularized quantile mixed model
#'
#' @keywords models
#' @param g the long format matrix of predictors (genetic factors) without intercept. Each row should be an observation vector.
#' @param y the vector of repeated measured responses. The current version of mixedBayes only supports continuous response.
#' @param e the long format matrix of environment (treatment) factors (a group of dummy variables).
#' @param X the long format matrix of the intercept and time effects (time effects are optional).
#' @param w the long format matrix of interactions between genetic factors and environment (treatment) factors.
#' @param k the number of repeated measurements.
#' @param iterations the number of MCMC iterations. The default value is 10,000.
#' @param burn.in the number of iterations for burn-in. If NULL, the first half of MCMC iterations will be used as burn-ins.
#' @param slope logical flag. If TRUE, random intercept-and-slope model will be used. Otherwise, random intercept model will be used. The default value is TRUE.
#' @param robust logical flag. If TRUE, robust methods will be used. Otherwise, non-robust methods will be used. The default value is TRUE.
#' @param quant the quantile level specified by users. The default value is 0.5.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be adopted to impose exact sparsity on regression coefficients. Otherwise, Laplacian shrinkage will be adopted. The default value is TRUE.
#' @param structure two choices are available. "bi-level" for selection on both the main and interaction effects corresponding to individual and group levels. "individual" for selections on individual-level only.
#' @return an object of class `mixedBayes' is returned, which is a list with component:
#' \item{posterior}{the posteriors of coefficients.}
#' \item{coefficient}{the estimated coefficients.}
#' \item{burn.in}{the total number of burn-ins.}
#' \item{iterations}{the total number of iterations.}
#'
#' @details Consider the data model described in "\code{\link{data}}":
#' \deqn{Y_{ij} = X_{ij}^\top\gamma_{0}+E_{ij}^\top\gamma_{1}+\sum_{l=1}^{p}G_{ijl}\gamma_{2l}+\sum_{l=1}^{p}W_{ijl}^\top\gamma_{3l}+Z_{ij}^\top\alpha_{i}+\epsilon_{ij}.}, with \eqn{W_{ij} = G_{ij}\bigotimes E_{ij}}.
#'
#' where \eqn{\gamma_{0}} is the coefficient vector for \eqn{X_{ij}}, \eqn{\gamma_{1}} is the coefficient vector for \eqn{E_{ij}}, \eqn{\gamma_{2l}} is the coefficient for the main effect of the \eqn{l}th genetic variant, and \eqn{\gamma_{3l}} is the coefficient vector for the interaction effect of the \eqn{l}th genetic variant with environment factors.
#'
#' For random intercept-and-slope model, \eqn{Z_{ij}^\top = (1,j)} and  \eqn{\alpha_{i} = (\alpha_{i1},\alpha_{i2})^\top}. For random intercept model, \eqn{Z_{ij}^\top = 1} and \eqn{\alpha_{i} = \alpha_{i1}}.
#'
#' When `structure="bi-level"`, bi-level selection will be conducted. If `structure="individual"`, individual-level selection will be conducted.
#'
#' When `slope=TRUE` (default), random intercept-and-slope model will be used as the mixed effects model.
#'
#' When `sparse=TRUE` (default), spike-and-slab priors are imposed to identify important main and interaction effects. Otherwise, Laplacian shrinkage will be used.
#'
#' When `robust=TRUE` (default), the distribution of \eqn{\epsilon_{ij}} is defined as an asymmetric Laplace distribution with density.
#'
#' \eqn{
#' f(\epsilon_{ij}|\theta,\tau) = \theta(1-\theta)\exp\left\{-\tau\rho_{\theta}(\epsilon_{ij})\right\}
#' }, (\eqn{i=1,\dots,n,j=1,\dots,k }), which leads to a Bayesian formulation of quantile regression. If `robust=FALSE`, \eqn{\epsilon_{ij}} follows a normal distribution.

#' \cr
#'
#' Please check the references for more details about the prior distributions.
#'
#' @seealso \code{\link{data}}
#'
#' @examples
#' data(data)
#'
#' ## default method (robust sparse bi-level selection under random intercept-and-slope model)

#' fit = mixedBayes(y,e,X,g,w,k,structure=c("bi-level"))
#' fit$coefficient
#'
#'## Compute TP and FP
#'b = selection(fit,sparse=TRUE)
#'index = which(coeff!=0)
#'pos = which(b != 0)
#'tp = length(intersect(index, pos))
#'fp = length(pos) - tp
#'list(tp=tp, fp=fp)

#' \donttest{
#' ## alternative: robust sparse individual level selections under random intercept-and-slope model
#' fit = mixedBayes(y,e,X,g,w,k,structure=c("individual"))
#' fit$coefficient
#'
#' ## alternative: non-robust sparse bi-level selection under random intercept-and-slope model
#' fit = mixedBayes(y,e,X,g,w,k,robust=FALSE, structure=c("bi-level"))
#' fit$coefficient
#'
#' ## alternative: robust sparse bi-level selection under random intercept model
#' fit = mixedBayes(y,e,X,g,w,k,slope=FALSE, structure=c("bi-level"))
#' fit$coefficient
#'
#' }
#'
#' @export

mixedBayes <- function(y,e,X,g,w,k, iterations=10000, burn.in=NULL, slope=TRUE, robust=TRUE, quant=0.5, sparse=TRUE, structure=c("bi-level","individual"))
{
  structure = match.arg(structure)

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
    out = LONRBGLSS(y,e,X,g,w,z,k,quant,iterations,sparse, structure)
  }else{
    out = LONBGLSS(y,e,X,g,w,z,k,iterations,sparse, structure)
  }
  }else{
    z = as.matrix(rep(1,k))
    if(robust){
      out = LONRBGLSS(y,e,X,g,w,z,k,quant,iterations,sparse, structure)
    }else{
      out = LONBGLSS(y,e,X,g,w,z,k,iterations,sparse, structure)
    }
  }

  coeff.treatment = apply(out$GS.gamma1[-(1:BI),,drop=FALSE], 2, stats::median);
  coeff.X = apply(out$GS.gamma0[-(1:BI),,drop=FALSE], 2, stats::median);
  coeff.main = apply(out$GS.gamma2[-(1:BI),,drop=FALSE], 2, stats::median);
  coeff.interaction = apply(out$GS.gamma3[-(1:BI),,drop=FALSE], 2, stats::median);
  coeff.random = apply(out$GS.alpha[-(1:BI),,drop=FALSE], 2, stats::median);
  coeff.random = matrix(coeff.random,ncol = length(y)/k)


  coefficient = list(treatment = coeff.treatment, X = coeff.X, main=coeff.main, interaction=coeff.interaction, random=coeff.random)

  fit = list(posterior = out, coefficient=coefficient, burn.in = BI, iterations=iterations)

  fit
}
