#' Probit analysis with random effects
#'
#' @description
#' Separate analysis over items.
#' @param formula Model formula, where multivariate responses may be given additively on the left hand side. Responses must be ordered factors.
#' @param subject Categorical variable encoding subjects.
#' @param data Date frame with data.
#' @param M Number of replications in simulation of E-step. Default: \code{M=10}.
#' @param dependence Text string (\code{"marginal"/"joint"}) deciding whether random effects are assumed independent or with a common joint normal distribution. Default: \code{dependence="marginal"}.
#' @note
#' A data frame must be provided, i.e. the \code{data} option is not optional. Variables not appearing in the \code{data}, but appearing in the \code{formula}, will be assumed to be random.
#' @return Object of class \code{probit}.
#' @export
probit <- function(formula,subject,data,M=10,dependence="marginal") {
  # work with tibbles
  data <- as_tibble(data)

  # find items and subjects
  items <- all.vars(formula[[2]])
  subjects <- unique(pull(data,!!enquo(subject)))

  # names of fixed and random effects
  vars <- all.vars(formula[[3]])
  fixed.eff <- intersect(vars,names(data))
  random.eff <- setdiff(vars,names(data))

  # add columns for random effects and initialize them at zero
  for (ii in random.eff) {
    data <- data %>% mutate(!! quo_name(enquo(ii)) := 0)
  }

  # data frames to contain lower and upper bounds
  alpha <- as.data.frame(matrix(NA,nrow(data),length(items)))
  beta  <- as.data.frame(matrix(NA,nrow(data),length(items)))
  names(alpha) <- items
  names(beta)  <- items

  # initial probit regressions with random effects fixated at zero
  regression <- vector("list",length(items))
  names(regression) <- items
  for (i in items) {
    if (!is.ordered(data[[i]])) stop(paste("Response variable",i,"must be an ordered factor"))
    regression[[i]] <- clm(eval(substitute(update(formula,y~.),list(y=as.name(i)))),
                           data=data,link="probit",na.action = na.exclude)
  }

  # Find alpha and beta from predictions
  # Hack: Set minimal and maximal value to -10 and 10, respectively.
  for (i in items) {
    tmp <- predict(regression[[i]],type="linear.predictor")
    alpha[[i]] <- pmax(-10,tmp$eta2)
    beta[[i]]  <- pmin(10,tmp$eta1)
  }

  # TO DO:
  for (i in subjects) {
    mydata <- data[pull(data,!!enquo(subject))==i,]
    subject.questions <- unique(pull(mydata,!!enquo(item)))
    for (ii in subject.questions) {

    }

  }

  # Find gamma from formula and parameter estimates

  # return
  return(list(fixed=fixed.eff,random=random.eff,res=res))
}
