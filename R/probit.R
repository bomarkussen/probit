#' Probit analysis with random effects
#'
#' @description
#' Separate analysis over item.
#' @param formula Model formula to be used separately for each item. Response variable \code{Y} must be a character stated as \code{ordered(S)} in the formula.
#' @param item Categorical variable encoding items.
#' @param subject Categorical variable encoding subjects.
#' @param data Date frame with data.
#' @param M Number of replications in simulation of E-step. Default: \code{M=10}.
#' @param dependence Text string (\code{"marginal"/"joint"}) deciding whether random effects are assumed independent or with a common joint normal distribution. Default: \code{dependence="marginal"}.
#' @note
#' A data frame must be provided, i.e. the \code{data} option is not optional. Variables not appearing in the \code{data}, but appearing in the \code{formula}, will be assumed to be random.
#' @return Object of class \code{probit}.
#' @export
probit <- function(formula,item,subject,data,M=10,dependence="marginal") {
  # work with tibbles
  data <- as_tibble(data)

  # only use complete cases
  if (sum(complete.cases(data))!=nrow(data)) {
    message(paste(nrow(data)-sum(complete.cases(data)),"rows with NA's has been removed from data."))
    data <- data[complete.cases(data),]
  }

  # make response variable into a character
  data[,all.vars(formula[[2]])] <- as.character(pull(data,all.vars(formula[[2]])))

  # find questions and subjects
  questions <- unique(pull(data,!!enquo(item)))
  subjects <- unique(pull(data,!!enquo(subject)))

  # names of fixed and random effects
  vars <- setdiff(all.vars(formula),as.character(formula[[2]]))
  fixed.eff <- intersect(vars,names(data))
  random.eff <- setdiff(vars,names(data))

  # add columns for random effects and initialize them at zero
  for (ii in random.eff) {
    data <- data %>% mutate(!! quo_name(enquo(ii)) := 0)
  }

  # add columns for lower and upper bounds, ie. alpha and beta
  if (any(is.element(names(data),c("lower_bound_alpha","upper_bound_beta")))) stop("Variable names 'lower_bound_alpha' and 'upper_bound_beta' are reserved, and may not appear in data.")
  data <- data %>% mutate(lower_bound_alpha=0,upper_bound_beta=0)

  # initial probit regressions with random effects fixated at zero
  res <- vector("list",length(questions))
  names(res) <- questions
  for (i in questions) {
    res[[i]] <- clm(formula,data=data[pull(data,!!enquo(item))==i,],link="probit")
  }

  # TO DO:
  # Find alpha and beta from predictions
  for (i in subjects) {
    mydata <- data[pull(data,!!enquo(subject))==i,]
    subject.questions <- unique(pull(mydata,!!enquo(item)))
    for (ii in subject.questions) {
      predict(res[[ii]],type="linear.predictor",newdata=mydata[pull(mydata,!!enquo(item))==ii,])

    }

  }

  # Find gamma from formula and parameter estimates

  # return
  return(list(fixed=fixed.eff,random=random.eff,res=res))
}
