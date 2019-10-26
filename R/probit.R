#' Probit analysis with random effects
#'
#' @description
#' Separate analysis over items.
#' @param formula Model formula, where multivariate responses may be given additively on the left hand side. Responses must be ordered factors.
#' @param subject Categorical variable encoding subjects.
#' @param data Date frame with data.
#' @param dependence Text string (\code{"marginal"/"joint"}) deciding whether random effects are assumed independent or with a common joint normal distribution. Default: \code{dependence="marginal"}.
#' @param Gamma Choleskey factor for initial variance of random effects. If \code{Gamma=NULL} then initialized at identity matrix. Default: \code{Gamma=NULL}.
#' @param M Number of replications in simulation of E-step. Default: \code{M=10}.
#' @param EMsteps Maximal number of EM steps. Default: \code{EMsteps=20}.
#' @param eps Absolute convergence criterion for total log likelihood. Default: \code{eps=1e-4}.
#' @param maxit Maximal number of steps in underlying coupling from the past algorithm. Default: \code{maxit=500}.
#' @param verbose If TRUE, then print convergence diagnostics. Default: \code{verbose=FALSE}.
#' @note
#' A data frame must be provided, i.e. the \code{data} option is not optional. Variables not appearing in the \code{data}, but appearing in the \code{formula}, will be assumed to be random.
#' @return Object of class \code{probit}.
#' @details
#' If \code{Gamma=NULL} then predicted random effects are initialized at zero. Otherwise they are simulated from normal distribution with variance \code{Gamma^T Gamma}.
#' @export
probit <- function(formula,subject,data,dependence="marginal",Gamma=NULL,M=10,EMsteps=20,eps=1e-4,maxit=500,verbose=FALSE) {
  # grab parameters ----

  # find items and subjects
  items <- all.vars(formula[[2]])
  subjects <- unique(data[[subject]])

  # names of fixed and random effects
  fixed.eff  <- intersect(all.vars(formula[[3]]),names(data))
  random.eff <-   setdiff(all.vars(formula[[3]]),names(data))

  # work with tibbles, make dataset as slim as possible, and make sanity check
  data <- as_tibble(data)
  data <- select(data,c(setdiff(all.vars(formula),random.eff),subject))
  if (is.element("my.clm.weight",names(data))) stop("'my.clm.weight' is used internally, and is not allowed as variable name")
  for (i in items) {
    if (!is.ordered(data[[i]])) stop(paste("Response variable",i,"must be an ordered factor"))
  }

  # zero.data = data with random effects fixated at zero
  # Remark: uses pred.random.eff as intermediate variable, which is reset below
  pred.random.eff <- matrix(0,length(subjects),1+length(random.eff))
  colnames(pred.random.eff) <- c(subject,random.eff)
  pred.random.eff <- as_tibble(pred.random.eff)
  pred.random.eff[,subject] <- subjects
  zero.data <- full_join(data,pred.random.eff,by=subject)


  # initialization ----

  # initialize prediction of random effects with M replicates
  pred.random.eff <- matrix(0,M*length(subjects),2+length(random.eff))
  colnames(pred.random.eff) <- c(subject,random.eff,"my.clm.weight")
  pred.random.eff <- as_tibble(pred.random.eff)
  pred.random.eff[,subject] <- rep(subjects,each=M)
  pred.random.eff[,"my.clm.weight"] <- 1/M

  # Cholesky factor for variance of random effects
  if (is.null(Gamma)) {
    # if not prespecified, then initialize as the identity matrix
    Gamma <- diag(nrow=length(random.eff))
  } else {
    # check format of Gamma
    if (nrow(Gamma)!=ncol(Gamma)) stop("Gamma must be NULL or a square matrix")
    if (nrow(Gamma)!=length(random.eff)) stop("Gamma must have number of columns equal to the number of random effects")
    # and then simulate
    pred.random.eff[,2:(1+length(random.eff))] <- as_tibble(matrix(rnorm(M*length(subjects)*length(random.eff)),M*length(subjects),length(random.eff))%*%Gamma)
  }

  # list to contain fit of probit regressions
  regression <- vector("list",length(items))
  names(regression) <- items

  # data frames to contain lower and upper bounds
  alpha <- as.data.frame(matrix(NA,nrow(data),length(items)))
  beta  <- as.data.frame(matrix(NA,nrow(data),length(items)))
  names(alpha) <- items
  names(beta)  <- items

  # array to contain slope against random effects
  gamma <- array(NA,dim=c(nrow(data),length(items),length(random.eff)))
  dimnames(gamma)[[2]] <- items
  dimnames(gamma)[[3]] <- random.eff

  # EM algorithm ----
  total.logLik <- -Inf
  for (iter in 0:EMsteps) {
    # M step ----

    # probit regressions at present values of random effects
    for (i in items) {
      regression[[i]] <- clm(eval(substitute(update(formula,y~.),list(y=as.name(i)))),
                             data=full_join(data,pred.random.eff,by=subject),
                             weights=my.clm.weight,
                             link="probit",na.action = na.exclude)
    }

    # variance matrix from predicted random effects
    # Remark: first after 0'th step
    if (iter > 0) {
      tmp <- as.matrix(pred.random.eff[,2:(1+length(random.eff))])
      if (dependence=="marginal") {
        Gamma <- diag(apply(tmp,2,function(x){sqrt(mean(x^2))}),nrow=length(random.eff))
      } else {
        Gamma <- chol(t(tmp)%*%tmp / nrow(pred.random.eff))
      }
    }

    # did likelihood improve
    tmp <- total.logLik
    total.logLik <- sum(sapply(regression,logLik))
    if (verbose) cat("EM step =",iter,": logLik =",total.logLik,"\n")
    if (total.logLik-eps < tmp) break

    # E step ----

    # Find alpha and beta from predictions
    # Remark: Set minimal and maximal value to -6 and 6, respectively.
    #         This corresponds to pnorm(-6) = 9.865876e-10. If we allow for
    #         smaller probabilities, then numerical instabilities appears to
    #         arise in the present implementation of the truncated simulations.
    for (i in items) {
      tmp <- predict(regression[[i]],type="linear.predictor",newdata = zero.data)
      alpha[[i]] <- pmax(-7,tmp$eta2)
      beta[[i]]  <- pmin(7,tmp$eta1)
    }

    # Find gamma via predictions
    # Remark: iii is introduced to avoid using the two extreme categories
    for (i in items) {
      p0 <- predict(regression[[i]],type="linear.predictor",newdata=zero.data)
      for (ii in random.eff) {
        iii <- (zero.data[[i]]==levels(data[[i]])[1])
        # predict with random effect increased by one
        # note: ordinal::clm() uses negative coefficients, and hence sign should not be changed
        zero.data[[ii]] <- 1
        p1 <- predict(regression[[i]],type="linear.predictor",newdata=zero.data)
        gamma[(!is.na(iii))&iii,i,ii]  <- p0$eta1[(!is.na(iii))&iii]  - p1$eta1[(!is.na(iii))&iii]
        gamma[(!is.na(iii))&!iii,i,ii] <- p0$eta2[(!is.na(iii))&!iii] - p1$eta2[(!is.na(iii))&!iii]
        # reset random effect at zero in order to restore zero.data
        zero.data[[ii]] <- 0
      }
    }

    # simulate predicted random effects
    for (i in subjects) {
      ii <- data[[subject]]==i
      tmp <- matrix(gamma[ii,,],sum(ii)*length(items),length(random.eff))
      iii <- apply(tmp,1,function(x){all(!is.na(x))})
      pred.random.eff[pred.random.eff[[subject]]==i,2:(1+length(random.eff))] <-
        r_truncated_multivariate_normal(n=M,
                                        alpha=c(as.matrix(alpha[ii,]))[iii],
                                        beta=c(as.matrix(beta[ii,]))[iii],
                                        gamma=t(tmp[iii,,drop=FALSE]),
                                        Gamma=Gamma,
                                        eps=1e-8,
                                        maxit=maxit)
    }

    # End of EM-loop
  }

  # return
  return(list(regression=regression,Gamma=Gamma))
}
