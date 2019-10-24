#' Probit analysis with random effects
#'
#' @description
#' Separate analysis over item.
#' @param formula Model formula to be used separately for each item.
#' @param item Categorical variable encoding items.
#' @param subject Categorical variable encoding subjects.
#' @param data Date frame with data.
#' @param dependence Text string (\code{"marginal"/"joint"}) deciding whether random effects are assumed independent or with a common joint normal distribution. Default: \code{dependence="marginal"}.
#' @note
#' A data frame must be provided, i.e. the \code{data} option is not optional. Variables not appearing in the \code{data}, but appearing in the \code{formula}, will be assumed to be random.
#' @return Object of class \code{probit}.
#' @export
probit <- function(formula,item,subject,data,dependence="marginal") {
  # find names of fixed and random effects
  vars <- setdiff(all.vars(formula),as.character(formula[[2]]))
  fixed.eff <- intersect(vars,names(data))
  random.eff <- setdiff(vars,names(data))

  # return
  return(list(fixed=fixed.eff,random=random.eff))
}
