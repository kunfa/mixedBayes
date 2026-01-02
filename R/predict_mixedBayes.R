#' Make predictions from a mixedBayes object
#'
#' @param object a mixedBayes object.
#' @param y a numeric vector of repeated-measure responses in long format. The observations
#'   should be stacked in the same order as the rows of X, e, g, and w.
#'   The current version only supports continuous response.
#' @param e the long-format design matrix for environment/treatment effects. In many applications,
#'   this is a set of dummy variables encoding treatment levels.
#' @param X the long-format design matrix for fixed effects, including an intercept and optionally
#'   time-related covariates.
#' @param g the long-format matrix of genetic predictors (main effects) without an intercept.
#' @param w the long-format matrix of gene-environment interaction terms. Typically, each row
#'   of w is constructed from the corresponding rows of g and e (e.g., via a Kronecker product),
#'   so that interaction terms represent products of genetic predictors with environment/treatment
#'   covariates.
#' @param k integer. Number of repeated measurements per subject.
#' @param slope logical flag. If TRUE, random intercept-and-slope model will be used.
#' @param loss character string specifying the prediction loss function.
#'        "L1" for mean absolute error;
#'        "L2" for mean squared error.
#'
#' @usage predict_mixedBayes(object, y, X, e, g, w, k, slope, loss)
#' @return an object of class `mixedBayes.pred' is returned, which is a list with components:
#' \item{pred_error}{prediction error.}
#' \item{y_hat}{predicted values of the repeated measured responses.}
#'
#' @rdname predict_mixedBayes
#' @seealso \code{\link{mixedBayes}}
#'
#' @examples
#' data(data)
#'
#' fit <- mixedBayes(y, e, X, g, w, k, structure = c("bi-level"))
#' pred1 <- predict_mixedBayes(fit, y, X, e, g, w, k, slope = TRUE, loss = "L1")
#' print(pred1$pred_error)
#' pred2 <- predict_mixedBayes(fit, y, X, e, g, w, k, slope = TRUE, loss = "L2")
#' print(pred2$pred_error)
#'
#' @export

predict_mixedBayes <- function(object, y, X, e, g, w, k,
                               slope = TRUE,
                               loss = c("L1", "L2")) {

  loss <- match.arg(loss)

  n  <- length(y)
  n1 <- n / k

  Zblock <- make_Zblock(k, n1, slope)

  y_hat <- X %*% object$coefficient$X +
    e %*% object$coefficient$treatment +
    g %*% object$coefficient$main +
    w %*% object$coefficient$interaction +
    Zblock %*% as.vector(object$coefficient$random)

  pred_error <- switch(
    loss,
    L1 = mean(abs(y - y_hat)),
    L2 = mean((y - y_hat)^2)
  )

  mixedBayes.pred <- list(
    y_hat = y_hat,
    pred_error = pred_error,
    loss = loss
  )

  class(mixedBayes.pred) <- "mixedBayes.pred"

  return(mixedBayes.pred)
}
