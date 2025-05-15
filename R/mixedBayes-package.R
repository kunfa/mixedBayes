#' @useDynLib mixedBayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

"_PACKAGE"
#' @keywords overview
#' @name mixedBayes-package
#' @title Bayesian Longitudinal Regularized Quantile Mixed Model
#' @aliases mixedBayes-package
#' @description In this package, we provide implementations of a set of high-dimensional robust Bayesian mixed-effect models to dissect longitudinal gene-environment interactions. The proposed method conducts robust Bayesian variable selection on both the main and interaction effects corresponding to individual and group levels (i.e. bi-level), respectively. Alternatively, selections only on individual levels by ignoring the grouping structure can also be performed. In addition, intra-cluster correlations among repeated measures are modeled via random intercept-and-slope and/or random intercept models. Imposing exact sparsity through spike-and-slab priors can be conducted on fixed effects with bi-level and/or individual level. In total, package mixedBayes provides implementations on 2 (robust and non-robust) × 2 ( types of fixed effects) × 2 ( types of random effects) × 2 (spike-and-slab or Laplacian priors) = 16 methods. Please read the details below for how to configure the method used.
#' @details The user friendly, integrated interface \strong{mixedBayes()} allows users to flexibly choose the fitting methods by specifying the following parameter:
#' \tabular{rl}{
#' slope: \tab whether to use random intercept-and-slope model or random intercept model.\cr\cr
#' robust: \tab whether to use robust or non-robust methods.\cr\cr
#' quant: \tab to specify different quantiles when using robust methods.\cr\cr
#' structure: \tab whether to specify bi-level or individual level.\cr\cr
#' sparse: \tab whether to use the spike-and-slab priors to impose sparsity.
#' }
#'
#' The function mixedBayes() returns a mixedBayes object that contains the posterior estimates of each coefficients.
#' S3 generic functions selection()and print() are implemented for mixedBayes objects.
#' selection() takes a mixedBayes object and returns the variable selection results.
#'
#' @references
#' Fan, K., Jiang, Y., Ma, S., Wang, W. and Wu, C. (2025). Robust Sparse Bayesian Regression for Longitudinal Gene-Environment Interactions.
#' {\emph{Journal of the Royal Statistical Society Series C: Applied Statistics}, qlaf027} \doi{10.1093/jrsssc/qlaf027}
#'
#' Zhou, F., Ren, J.,  Li, G., Jiang, Y., Li, X., Wang, W. and Wu, C. (2019). Penalized Variable Selection for Lipid-Environment Interactions in a Longitudinal Lipidomics Study.
#' {\emph{Genes}, 10(12), 1002} \doi{10.3390/genes10121002}
#'
#' Zhou, F., Ren, J., Liu, Y., Li, X., Wang, W., and Wu, C. (2022). Interep: An r package for high-dimensional interaction analysis of the repeated measurement data.
#' {\emph{Genes}, 13(3), 544} \doi{10.3390/genes13030544}
#'
#' Zhou, F., Lu, X., Ren, J., Fan, K., Ma, S., and Wu, C. (2022). Sparse group variable selection for gene–environment interactions in the longitudinal study.
#' {\emph{Genetic epidemiology}, 46(5-6), 317-340} \doi{10.1002/gepi.22461}
#'
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2023). Robust Bayesian variable selection for gene-environment interactions.
#' {\emph{Biometrics},79(2),684-694 } \doi{10.1111/biom.13670}
#'
#' Wu, C., and Ma, S. (2015). A selective review of robust variable selection with applications in bioinformatics.
#' {\emph{Briefings in Bioinformatics}, 16(5), 873–883} \doi{10.1093/bib/bbu046}
#'
#' Zhou, F., Ren, J., Lu, X., Ma, S. and Wu, C. (2021). Gene–Environment Interaction: a Variable Selection Perspective.
#' {\emph{Epistasis. Methods in Molecular Biology.} 2212:191–223} \doi{10.1007/978-1-0716-0947-7_13}
#'
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y. and Wu, C. (2020) Semi-parametric Bayesian variable selection for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 39: 617– 638} \doi{10.1002/sim.8434}
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y. and Ma, S. (2018). Dissecting gene-environment interactions: A penalized robust approach accounting for hierarchical structures.
#' {\emph{Statistics in Medicine}, 37:437–456} \doi{10.1002/sim.7518}
#'
#' Wu, C., Cui, Y., and Ma, S. (2014). Integrative analysis of gene–environment interactions under a multi–response partially linear varying coefficient model.
#' {\emph{Statistics in Medicine}, 33(28), 4988–4998} \doi{10.1002/sim.6287}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2013). High dimensional variable selection for gene-environment interactions.
#' {\emph{Technical Report. Michigan State University.}}
#'
#' @seealso \code{\link{mixedBayes}}
NULL
