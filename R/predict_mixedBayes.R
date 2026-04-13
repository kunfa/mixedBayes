#' Make predictions from a mixedBayes object
#'
#' @param object a mixedBayes object.
#' @param y a numeric vector of repeated-measure responses in long format.
#'   The current version only supports continuous response.
#' @param e the long-format design matrix for environment/treatment effects. In applications,
#'   this is a group of dummy variables encoding treatment levels.
#' @param X the long-format design matrix, including an intercept and optionally
#'   time-related covariates.
#' @param g the long-format matrix of genetic predictors.
#' @param k integer. Number of repeated measurements per subject.
#' @param slope logical flag. If TRUE, random intercept-and-slope model will be used.
#' @param loss character string specifying the prediction loss function.
#'        "L1" for mean absolute error;
#'        "L2" for mean squared error.
#'
#' @usage predict_mixedBayes(object, y, X, e, g, k, slope, loss)
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
#' fit <- mixedBayes(y, e, X, g, k, structure = "bilevel")
#' pred1 <- predict_mixedBayes(fit, y, X, e, g, k, slope = TRUE, loss = "L1")
#' print(pred1$pred_error)
#' fit <- mixedBayes(y, e, X, g, k, robust =FALSE, quant =NULL,structure = "bilevel")
#' pred2 <- predict_mixedBayes(fit, y, X, e, g, k, slope = TRUE, loss = "L2")
#' print(pred2$pred_error)
#'
#' @export

predict_mixedBayes <- function(object, y, X, e, g, k,
                               slope = TRUE,
                               loss = "L1") {
  if (length(k) != 1 || !is.numeric(k) || is.na(k) || k <= 0 || k %% 1 != 0) {
    stop("k must be a positive integer.")
  }
  if (!is.logical(slope) || length(slope) != 1 || is.na(slope)) {
    stop("slope must be TRUE or FALSE.")
  }

  n  <- length(y)
  if (n %% k != 0) {
    stop("length(y) must be divisible by k.")
  }

  if (nrow(X) != n || nrow(e) != n || nrow(g) != n) {
    stop("X, e, and g must all have the same number of rows as length(y).")
  }

  n1 <- n / k

  Zblock <- make_Zblock(k, n1, slope)
  m = ncol(g)
  w = c()

  for (i in 1:m)
  {

    w = cbind(w,g[,i]*e)

  }

  y_hat <- X %*% object$coefficient$X +
    e %*% object$coefficient$treatment +
    g %*% object$coefficient$main +
    w %*% object$coefficient$interaction +
    Zblock %*% as.vector(object$coefficient$random)

  if (loss == "L1") {
    pred_error <- mean(abs(y - y_hat))
  } else if (loss == "L2") {
    pred_error <- mean((y - y_hat)^2)
  } else {
    stop("loss must be either 'L1' or 'L2'.")
  }

  mixedBayes.pred <- list(
    y_hat = y_hat,
    pred_error = pred_error,
    loss = loss
  )

  class(mixedBayes.pred) <- "mixedBayes.pred"

  return(mixedBayes.pred)
}
