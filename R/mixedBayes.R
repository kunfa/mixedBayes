#' fit a Bayesian longitudinal regularized quantile mixed model
#'
#' @keywords models
#' @param y a numeric vector of repeated-measure responses in long format.
#'   The current version only supports continuous response.
#' @param e the long-format design matrix for environment/treatment effects. In applications,
#'   this is a set of dummy variables encoding treatment levels.
#' @param X the long-format design matrix, including an intercept and optionally
#'   time-related covariates.
#' @param g the long-format matrix of genetic predictors.
#' @param w the long-format matrix of gene-environment interaction terms.
#' @param k integer. Number of repeated measurements per subject.
#' @param iterations the number of MCMC iterations. The default value is 10,000.
#' @param burn.in the number of iterations for burn-in. If NULL, no burn-in is applied and all MCMC samples are retained.
#' @param slope logical flag. If TRUE, random intercept-and-slope model will be used. Otherwise, random intercept model will be used. The default value is TRUE.
#' @param robust logical flag. If TRUE, robust methods will be used. Otherwise, non-robust methods will be used. The default value is TRUE.
#' @param quant the quantile level specified by users. The default value is 0.5.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be adopted to impose exact sparsity on regression coefficients. Otherwise, Laplacian shrinkage will be adopted. The default value is TRUE.
#' @param structure two choices are available. "bi-level" for selection on both the main and interaction effects corresponding to individual and group levels. "individual" for selections on individual-level only.
#' @return an object of class `mixedBayes' is returned, which is a list with component:
#' \item{posterior}{posterior samples for fixed effects and random effects.}
#' \item{coefficient}{posterior median estimates of coefficients for fixed effects and random effects.}
#'
#' @details
#' \strong{Data layout}
#'
#' Consider a longitudinal study with repeated measurements per subject. The response vector y
#' and the design matrices X, e, g, and w must all be provided in long format and share the same
#' row ordering. In practice, each row corresponds to one observation from a particular subject at a
#' particular time point.
#'
#' \strong{Model}
#'
#' Consider the data model described in "\code{\link{data}}":
#' \deqn{Y_{ij} = X_{ij}^\top\gamma_{0}+E_{ij}^\top\gamma_{1}+\sum_{l=1}^{p}G_{ijl}\gamma_{2l}+\sum_{l=1}^{p}W_{ijl}^\top\gamma_{3l}+Z_{ij}^\top\alpha_{i}+\epsilon_{ij}.}, with \eqn{W_{ij} = G_{ij}\bigotimes E_{ij}}.
#'
#' where \eqn{\gamma_{0}} is the coefficient vector for \eqn{X_{ij}}, \eqn{\gamma_{1}} is the coefficient vector for \eqn{E_{ij}}, \eqn{\gamma_{2l}} is the coefficient for the main effect of the \eqn{l}th genetic variant, and \eqn{\gamma_{3l}} is the coefficient vector for the interaction effect of the \eqn{l}th genetic variant with environment factors.
#'
#' The subject-specific random effects \eqn{\alpha_i} capture within-subject correlation. For random intercept-and-slope model, \eqn{Z_{ij}^\top = (1,j)} and  \eqn{\alpha_{i} = (\alpha_{i1},\alpha_{i2})^\top}. For random intercept model, \eqn{Z_{ij}^\top = 1} and \eqn{\alpha_{i} = \alpha_{i1}}.
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
  X = as.matrix(X)
  structure = match.arg(structure)

  if(slope){
    k_1 = c(1:k)
    k_1 = k_1 - mean(k_1)
    z = cbind(rep(1,k),k_1)
  if(robust){
    out = LONRBGLSS(y,e,X,g,w,z,k,quant,iterations,sparse, structure, iterations,burn.in)
  }else{
    out = LONBGLSS(y,e,X,g,w,z,k,iterations,sparse, structure, iterations,burn.in)
  }
  }else{
    z = as.matrix(rep(1,k))
    if(robust){
      out = LONRBGLSS_1(y,e,X,g,w,z,k,quant,iterations,sparse, structure, iterations,burn.in)
    }else{
      out = LONBGLSS_1(y,e,X,g,w,z,k,iterations,sparse, structure, iterations,burn.in)
    }
  }

  coeff.treatment = apply(out$GS.gamma1, 2, stats::median);
  coeff.X = apply(as.matrix(out$GS.gamma0), 2, stats::median);
  coeff.main = apply(out$GS.gamma2, 2, stats::median);
  coeff.interaction = apply(out$GS.gamma3, 2, stats::median);
  coeff.random = apply(out$GS.alpha, 2, stats::median);
  coeff.random = matrix(coeff.random,ncol = length(y)/k)


  coefficient = list(treatment = coeff.treatment, X = coeff.X, main=coeff.main, interaction=coeff.interaction, random=coeff.random)

  fit = list(posterior = out, coefficient=coefficient)

  fit
}
