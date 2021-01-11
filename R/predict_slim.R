#' Prediction from new dataset
#'
#' @description
#' Predictions from a linear model.
#' @param object of clasee \code{\link[stats]{lm}} or \code{\link[biglm]{biglm}}.
#' @param newdata data frame on which predictions should be made.
#'
#' @return vector with predictions.
#'
#' @details
#' A slim implementation of predictions from a linear model
#' on a new dataset. This is made in replacement of \code{\link[stats]{predict}}
#' and \code{\link[biglm]{predict.bigglm}}, which do not work on
#' \code{\link[biglm]{biglm}}-objects or does not include a possible
#' \code{\link[stats]{offset}}, respectively.
#'
#' @note The implementation uses \code{biglm:::coef.biglm} as this function
#' is not exported from the \code{\link[biglm]{biglm}}-package.
#'
#' @export
predict_slim <- function(object,newdata) {
  if (class(object)=="lm") {theta <- coef(object)} else {theta <- biglm:::coef.biglm(object)}

  theta_notNA <- !is.na(theta)
  X <- model.matrix(delete.response(terms(object)),
                    model.frame(delete.response(terms(object)),newdata))
  res <- c(X[,theta_notNA[colnames(X)],drop=FALSE] %*%
           theta[is.element(names(theta),colnames(X)) & theta_notNA])

  if (!is.null(attr(object$terms, "offset"))) res <- res +
      model.offset(model.frame(delete.response(terms(object)),newdata))

  # return
  return(res)
}
