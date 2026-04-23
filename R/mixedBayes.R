#' fit a Bayesian longitudinal regularized quantile mixed model
#'
#' @keywords models
#' @param y a numeric vector of repeated-measure responses in long format. The current version only supports continuous response.
#' @param e the long-format matrix for environmental/treatment effects. In applications, this corresponds to a set of dummy variables encoding treatment levels.
#' @param X the long-format matrix, including an intercept term and optionally time covariates.
#' @param g the long-format matrix of genetic predictors.
#' @param k a positive integer. Number of repeated measurements per subject.
#' @param iterations the number of MCMC iterations. The default value is 10,000.
#' @param burn.in the number of iterations for burn-in. If NULL, no burn-in is applied and all MCMC samples are retained. The default value is 5,000.
#' @param slope logical flag. If TRUE, random intercept-and-slope model will be used. Otherwise, random intercept model will be used. The default value is TRUE.
#' @param robust logical flag. If TRUE, robust methods will be used. Otherwise, non-robust methods will be used. The default value is TRUE.
#' @param quant the quantile level specified by users. Required when robust = TRUE. Ignored (set to NULL) when robust = FALSE. The default value is 0.5.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be adopted to impose exact sparsity on regression coefficients. Otherwise, Laplacian shrinkage will be adopted. The default value is TRUE.
#' @param structure two choices are available. "bilevel" performs selection on both main effects and interaction effects corresponding to individual and group levels, whereas "individual" performs selections only on individual levels by ignoring the group structure.
#' @return An object of class \code{"mixedBayes"} is returned, which is a list with components:
#' \item{posterior}{Posterior samples for fixed effects and random effects.}
#' \item{coefficient}{Posterior median estimates of fixed effects and random effects.}
#' \item{sparse}{Logical value indicating whether spike-and-slab priors were used.}
#' \item{robust}{Logical value indicating whether the robust (quantile regression) model was used.}
#' \item{slope}{Logical value indicating whether a random intercept-and-slope model was used.}
#' \item{quant}{Quantile level used when \code{robust = TRUE}; otherwise \code{NULL}.}
#' \item{structure}{Character string specifying the selection structure ("bilevel" or "individual").}
#' \item{k}{Number of repeated measurements per subject.}
#' \item{iterations}{Total number of MCMC iterations.}
#' \item{burn.in}{Number of burn-in iterations.}
#'
#' @details
#' \strong{Data layout}
#'
#' Consider a longitudinal study with repeated measurements per subject, where each subject has the same number of repeated measurements.
#'
#' The response vector y and the design matrices X, e and g must all be provided in long format and share the same
#' row ordering. In practice, each row corresponds to one observation from a particular subject at a
#' particular time point.
#'
#' The interaction terms between genetic and environmental factors (G×E) are constructed internally
#' within the function and therefore do not need to be provided by the user.
#'
#'
#' \strong{Model}
#'
#' Consider the data model described in "\code{\link{data}}":
#' \deqn{Y_{ij} = \boldsymbol{X^\top_{ij}}\boldsymbol{\gamma_{0}}+\boldsymbol{E^\top_{ij}}\boldsymbol{\gamma_{1}}+\boldsymbol{G^\top_{ij}}\boldsymbol{\gamma_{2}}+(\boldsymbol{G_{ij}}\bigotimes \boldsymbol{E_{ij}})^\top\boldsymbol{\gamma_{3}}+\boldsymbol{Z^\top_{ij}}\boldsymbol{\alpha_{i}}+\epsilon_{ij}.}
#'
#' Here \eqn{\boldsymbol{\gamma_{0}}} is the coefficient vector for \eqn{\boldsymbol{X_{ij}}}, \eqn{\boldsymbol{\gamma_{1}}} is the coefficient vector for \eqn{\boldsymbol{E_{ij}}}, \eqn{\boldsymbol{\gamma_{2}}} is the coefficient vector for the genetic variants, and \eqn{\boldsymbol{\gamma_{3}}} is the coefficient vector for the interactions of the genetic variants with environment factors.
#'
#' where \eqn{\boldsymbol{\gamma_{1}}=(\gamma_{11},\dots,\gamma_{1p})^\top}, \eqn{\boldsymbol{\gamma_{2}}=(\gamma_{21},\dots,\gamma_{2m})^\top}, \eqn{\boldsymbol{\gamma_{3}}=(\boldsymbol{\gamma_{31}},\dots,\boldsymbol{\gamma_{3m}})^\top} where \eqn{\boldsymbol{\gamma_{3l}}=(\gamma_{3l1},\dots,\gamma_{3lp})^\top} for \eqn{l=1,\dots,m}.
#'
#' The subject-specific random effects \eqn{\boldsymbol{\alpha_i}} capture within-subject correlation. For random intercept-and-slope model, \eqn{\boldsymbol{Z^\top_{ij}} = (1,j)} and  \eqn{\boldsymbol{\alpha_{i}} = (\alpha_{i1},\alpha_{i2})^\top}. For random intercept model, \eqn{Z_{ij} = 1} and \eqn{\alpha_{i} = \alpha_{i1}}.
#'
#' When `structure=bilevel`(default), bi-level selection on main and interaction effects will be conducted corresponding to individual and group levels, respectively. When `structure="individual"`, selections only on individual levels by ignoring the group structure will be performed.
#'
#' When `slope=TRUE` (default), random intercept-and-slope model will be used as the mixed effects model. Otherwise, random intercept model will be used.
#'
#' When `sparse=TRUE` (default), spike-and-slab priors are imposed to identify important main and interaction effects. Otherwise, Laplacian shrinkage will be used.
#'
#' When `robust=TRUE` (default), the distribution of \eqn{\epsilon_{ij}} is defined as an asymmetric Laplace distribution with density.
#'
#' \eqn{
#' f(\epsilon_{ij}|\theta,\tau) = \theta(1-\theta)\exp\left\{-\tau\rho_{\theta}(\epsilon_{ij})\right\}
#' }, (\eqn{i=1,\dots,n,j=1,\dots,k }), which leads to a Bayesian formulation of quantile regression. Otherwise, \eqn{\epsilon_{ij}} follows a normal distribution.

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

#' fit = mixedBayes(y,e,X,g,k,structure="bilevel")
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
#' fit = mixedBayes(y,e,X,g,k,structure="individual")
#' fit$coefficient
#'
#' ## alternative: non-robust sparse bi-level selection under random intercept-and-slope model
#' fit = mixedBayes(y,e,X,g,k,robust=FALSE, quant = NULL, structure="bilevel")
#' fit$coefficient
#'
#' ## alternative: robust sparse bi-level selection under random intercept model
#' fit = mixedBayes(y,e,X,g,k,slope=FALSE, structure="bilevel")
#' fit$coefficient
#'
#' }
#'
#' @export

mixedBayes <- function(y, e, X, g, k,
                       iterations = 10000,
                       burn.in = 5000,
                       slope = TRUE,
                       robust = TRUE,
                       quant = 0.5,
                       sparse = TRUE,
                       structure = "bilevel") {

  X <- as.matrix(X)
  e <- as.matrix(e)
  g <- as.matrix(g)
  y <- as.numeric(y)

  n <- length(y)
  E <- cbind(e, X)
  o <- ncol(X)
  q <- ncol(E)

  if (nrow(X) != n || nrow(e) != n || nrow(g) != n) {
    stop("X, e, and g must all have the same number of rows as length(y).")
  }

  if (length(k) != 1 || !is.numeric(k) || is.na(k) || k <= 0 || k %% 1 != 0) {
    stop("k must be a positive integer.")
  }

  if (n %% k != 0) {
    stop("length(y) must be divisible by k.")
  }

  if (!is.logical(slope) || length(slope) != 1 || is.na(slope)) {
    stop("slope must be TRUE or FALSE.")
  }

  if (!is.logical(robust) || length(robust) != 1 || is.na(robust)) {
    stop("robust must be TRUE or FALSE.")
  }

  if (!is.logical(sparse) || length(sparse) != 1 || is.na(sparse)) {
    stop("sparse must be TRUE or FALSE.")
  }

  if (length(structure) != 1 || !is.character(structure) || is.na(structure) ||
      !(structure %in% c("bilevel", "individual"))) {
    stop("structure must be either 'bilevel' or 'individual'.")
  }

  if (length(iterations) != 1 || !is.numeric(iterations) || is.na(iterations) ||
      iterations <= 0 || iterations %% 1 != 0) {
    stop("iterations must be a positive integer.")
  }

  if (is.null(burn.in)) {
    BI <- 0L
  } else if (length(burn.in) == 1 && is.numeric(burn.in) && !is.na(burn.in) &&
             burn.in >= 0 && burn.in %% 1 == 0) {
    BI <- as.integer(burn.in)
  } else {
    stop("burn.in must be NULL or a nonnegative integer.")
  }

  if (iterations <= BI) {
    stop("iterations must be larger than burn.in.")
  }

  if (robust) {
    if (is.null(quant) || length(quant) != 1 || !is.numeric(quant) || is.na(quant) ||
        quant <= 0 || quant >= 1) {
      stop("quant must be a single number between 0 and 1 when robust = TRUE.")
    }
  } else {
    if (!is.null(quant)) {
      stop("quant must be NULL when robust = FALSE.")
    }
  }

  if (slope) {
    k_1 <- seq_len(k)
    k_1 <- k_1 - mean(k_1)
    z <- cbind(rep(1, k), k_1)

    if (robust) {
      fit0 <- LONRBGLSS(y, e, X, g, z, k, quant, iterations, sparse, structure, iterations, burn.in)
    } else {
      fit0 <- LONBGLSS(y, e, X, g, z, k, iterations, sparse, structure, iterations, burn.in)
    }

  } else {
    z <- matrix(1, nrow = k, ncol = 1)

    if (robust) {
      fit0 <- LONRBGLSS_1(y, e, X, g, z, k, quant, iterations, sparse, structure, iterations, burn.in)
    } else {
      fit0 <- LONBGLSS_1(y, e, X, g, z, k, iterations, sparse, structure, iterations, burn.in)
    }
  }

  if (BI == 0) {
    out <- list(
      GS.gamma1 = fit0$GS.alpha[, 1:(q - o), drop = FALSE],
      GS.gamma0 = fit0$GS.alpha[, -(1:(q - o)), drop = FALSE],
      GS.gamma2 = fit0$GS.beta,
      GS.gamma3 = fit0$GS.eta,
      GS.alpha  = fit0$GS.ata
    )
  } else {
    burn_rows <- seq_len(BI)
    out <- list(
      GS.gamma1 = fit0$GS.alpha[-burn_rows, 1:(q - o), drop = FALSE],
      GS.gamma0 = fit0$GS.alpha[-burn_rows, -(1:(q - o)), drop = FALSE],
      GS.gamma2 = fit0$GS.beta[-burn_rows, , drop = FALSE],
      GS.gamma3 = fit0$GS.eta[-burn_rows, , drop = FALSE],
      GS.alpha  = fit0$GS.ata[-burn_rows, , drop = FALSE]
    )
  }

  coeff.treatment   <- apply(out$GS.gamma1, 2, stats::median)
  coeff.X           <- apply(out$GS.gamma0, 2, stats::median)
  coeff.main        <- apply(out$GS.gamma2, 2, stats::median)
  coeff.interaction <- apply(out$GS.gamma3, 2, stats::median)
  coeff.random      <- apply(out$GS.alpha, 2, stats::median)
  coeff.random      <- matrix(coeff.random, ncol = n / k)

  coefficient <- list(
    treatment   = coeff.treatment,
    X           = coeff.X,
    main        = coeff.main,
    interaction = coeff.interaction,
    random      = coeff.random
  )

  fit <- list(
    posterior   = out,
    coefficient = coefficient,
    sparse      = sparse,
    robust      = robust,
    slope       = slope,
    quant       = quant,
    structure   = structure,
    k           = k,
    iterations  = iterations,
    burn.in     = BI
  )

  class(fit) <- "mixedBayes"

  return(fit)
}
