#' simulated data for demonstrating the features of rolong
#'
#' Simulated gene expression data for demonstrating the features of rolong.
#'
#' @docType data
#' @keywords datasets
#' @name data
#' @aliases data y e g w k
#' @format The data object consists of five components: y, e, g, w and k.
#'
#' @details
#'
#' \strong{The data and model setting}
#'
#' Consider a longitudinal study on \eqn{n} subjects with \eqn{k} repeated measurement for each subject. Let \eqn{Y_{ij}} be the measurement for the \eqn{i}th subject at each time point \eqn{j}(\eqn{1\leq i \leq n, 1\leq j \leq k}) .We use a \eqn{m}-dimensional vector \eqn{G_{ij}} to denote the genetics factors, where\eqn{L_{ij} = (L_{ij1},...,L_{ijm})^{T}}. Also, we use \eqn{p}-dimensional vector \eqn{E_{ij}} to denote the treatment factors, where \eqn{T_{ij} = (T_{ij1},...,T_{ijp})^{T}}. \eqn{X_{ij} = (1, j, j^2)^{T}}. \eqn{Z_{ij}} is a \eqn{h \times 1} covariate associated with random effects and \eqn{\alpha_{i}} is a \eqn{h\times 1} vector of random effects.  At the beginning,  the interaction effects is modeled as the product of genomics features and treatment factor with 4 different levels. After representing the treatment factors as three dummy variables, the identification of the lipid by treatment interaction needs to be performed as group level.  Combing the genetics factors, treatment factors and their interactions that associated with the longitudinal phenotype, we have the following mixed-effects model:
#' \deqn{Y_{ij} = X_{ij}^{T}\gamma_{0}+E_{ij}^{T}\gamma_{1}+G_{ij}^{T}\gamma_{2}+(G_{ij}\bigotimes E_{ij})^{T}\gamma_{3}+Z_{ij}^{T}\alpha_{i}+\epsilon_{ij}.}
#' where \eqn{\gamma_{1}},\eqn{\gamma_{2}},\eqn{\gamma_{3}} are \eqn{p},\eqn{m} and \eqn{mp} dimensional vectors that represent the coefficients of the treatment effects, the genetics effects and interactions effects, respectively. Accommodating the Kronecker product of the \eqn{m} - dimensional vector \eqn{L_{ij}} and the \eqn{p}-dimensional vector \eqn{E_{ij}}, the interactions between genetics and treatment factors can be expressed as a \eqn{mp}-dimensional vector, denoted as the following form:
#' \deqn{G_{ij}\bigotimes E_{ij} = [E_{ij1}E_{ij1},E_{ij2}E_{ij2},...,E_{ij1}E_{ijp},E_{ij2}E_{ij1},...,E_{ijm}E_{ijp}]^{T}.}
#' When \eqn{h=1}, the model becomes a mixed-effects model with random intercept only.
#'
#' @examples
#' data(data)
#' dim(y)
#' dim(g)
#' dim(e)
#' dim(w)
#' print(k)
#' print(coeff)
#'
#' @seealso \code{\link{rolong}}
NULL
