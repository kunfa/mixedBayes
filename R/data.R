#' simulated data for demonstrating the features of mixedBayes
#'
#' Simulated gene expression data for demonstrating the features of mixedBayes.
#'
#' @docType data
#' @keywords datasets
#' @name data
#' @aliases data y e X g w k coeff
#' @format The data object consists of seven components: y, e, X, g, w ,k and coeff. coeff contains the true values of parameters (main and interaction effects) used for generating Y.
#'
#' @details
#'
#' \strong{The data and model setting}
#'
#' Consider a longitudinal study on \eqn{n} subjects with \eqn{k} repeated measurement for each subject. Let \eqn{Y_{ij}} be the measurement for the \eqn{i}th subject at each time point \eqn{j}(\eqn{1\leq i \leq n, 1\leq j \leq k}) .We use the \eqn{m}-dimensional vector \eqn{G_{ij}} to denote measurements of genetics factors for the \eqn{i}th subject at time point \eqn{j}, where \eqn{G_{ij} = (G_{ij1},...,G_{ijm})^\top}. Also, we use \eqn{p}-dimensional vector \eqn{E_{ij}} to denote the environment factors, where \eqn{E_{ij} = (E_{ij1},...,E_{ijp})^\top}. \eqn{X_{ij} = (1, T_{ij})^\top}, where \eqn{T_{ij}^\top} is a vector of time effects . \eqn{Z_{ij}} is a \eqn{h \times 1} covariate associated with random effects and \eqn{\alpha_{i}} is a \eqn{h\times 1} vector of random effects. In a typical one-way repeated measure ANOVA with a fixed number (say four) of factor levels, the environment (or treatment) factor is modelled as a group of three dummy variables. Therefore, gene-environment (or treatment) interaction leads to variable selections on individual levels (main effects) and group levels (interaction effect) simultaneously. Considering the genetics factors, environment (or treatment) factors and their interactions that are jointly associated with the longitudinal phenotype, we have the following mixed-effects model:
#' \deqn{Y_{ij} = X_{ij}^\top\gamma_{0}+E_{ij}^\top\gamma_{1}+G_{ij}^\top\gamma_{2}+(G_{ij}\bigotimes E_{ij})^\top\gamma_{3}+Z_{ij}^\top\alpha_{i}+\epsilon_{ij}.}
#' where \eqn{\gamma_{1}},\eqn{\gamma_{2}},\eqn{\gamma_{3}} are \eqn{p},\eqn{m} and \eqn{mp} dimensional vectors that represent the coefficients of the environment effects, the genetics effects and interactions effects, respectively. In addition, \eqn{\gamma_0} is the coefficient vector for \eqn{X_{ij}}.
#' The gene–environment interactions that can be expressed as a Kronecker product between the two types of main effects as a \eqn{mp}-dimensional vector:
#' \deqn{G_{ij}\bigotimes E_{ij} = [G_{ij1}E_{ij1},G_{ij1}E_{ij2},...,G_{ij1}E_{ijp},G_{ij2}E_{ij1},...,G_{ijm}E_{ijp}]^\top.}
#'
#' The above model also includes \eqn{Z_{ij}} with random effects \eqn{\alpha_{i}} to account for intra-correlations among repeated measurements.
#' For random intercept-and-slope model, \eqn{Z_{ij}^\top = (1,j)} and  \eqn{\alpha_{i} = (\alpha_{i1},\alpha_{i2})^\top}. For random intercept model, \eqn{Z_{ij}^\top = 1} and \eqn{\alpha_{i} = \alpha_{i1}}.
#'
#' In this example, we generate data under random intercept-and-slope model.
#'
#' @examples
#' data(data)
#' length(y)
#' dim(g)
#' dim(e)
#' dim(w)
#' print(k)
#' print(X)
#' print(coeff)
#'
#' @seealso \code{\link{mixedBayes}}
NULL
