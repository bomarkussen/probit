#' @title Simulate random effects
#'
#' @description
#' Simulate random effects from a probit object.
#' @param x Object of class \code{\link{probit-class}}
#' @param maxit Maximal number of steps in coupling from the past algorithm. Default: \code{maxit=500}.
#' @export
ranef.probit <- function(x,maxit=500) {
  # find subjects
  subjects <- unique(x$data[[x$subject]])

  # zero.data = data with random effects fixated at zero
  zero.data <- matrix(0,length(subjects),1+length(x$random.eff))
  colnames(zero.data) <- c(x$subject,x$random.eff)
  zero.data <- as_tibble(zero.data)
  zero.data[,x$subject] <- subjects
  zero.data <- full_join(x$data,zero.data,by=x$subject)

  # Find alpha and beta from predictions
  # Remark: Set minimal and maximal value to -6 and 6, respectively.
  #         This corresponds to pnorm(-6) = 9.865876e-10.
  alpha <- as_tibble(matrix(NA,nrow(x$data),length(x$items)))
  beta  <- as_tibble(matrix(NA,nrow(x$data),length(x$items)))
  names(alpha) <- x$items
  names(beta)  <- x$items
  for (i in x$items) {
    tmp <- predict(x$regression[[i]],type="linear.predictor",newdata = zero.data)
    alpha[[i]] <- pmax(-6,tmp$eta2)
    beta[[i]]  <- pmin(6,tmp$eta1)
  }

  # Find gamma via predictions
  # Remark: iii is introduced to avoid using the two extreme categories
  gamma <- array(NA,dim=c(nrow(x$data),length(x$items),length(x$random.eff)))
  dimnames(gamma)[[2]] <- x$items
  dimnames(gamma)[[3]] <- x$random.eff
  for (i in x$items) {
    p0 <- predict(x$regression[[i]],type="linear.predictor",newdata=zero.data)
    for (ii in x$random.eff) {
      iii <- (zero.data[[i]]==levels(x$data[[i]])[1])
      # predict with random effect increased by one
      zero.data[[ii]] <- 1
      p1 <- predict(x$regression[[i]],type="linear.predictor",newdata=zero.data)
      gamma[(!is.na(iii))&iii,i,ii]  <- p0$eta1[(!is.na(iii))&iii]  - p1$eta1[(!is.na(iii))&iii]
      gamma[(!is.na(iii))&!iii,i,ii] <- p0$eta2[(!is.na(iii))&!iii] - p1$eta2[(!is.na(iii))&!iii]
      # reset random effect at zero in order to restore zero.data
      zero.data[[ii]] <- 0
    }
  }

  # set-up tibble to contain residuals
  residual <- matrix(NA,nrow(x$data),length(x$items))
  colnames(residual) <- x$items
  residual <- as_tibble(residual)

  # set-up tibble to contain random effects
  Z <- matrix(0,length(subjects),1+length(x$random.eff))
  colnames(Z) <- c(x$subject,x$random.eff)
  Z <- as_tibble(Z)
  Z[,x$subject] <- subjects

  # simulate from estimated models
  for (i in subjects) {
    ii <- x$data[[x$subject]]==i
    tmp <- matrix(gamma[ii,,],sum(ii)*length(x$items),length(x$random.eff))
    iii <- apply(tmp,1,function(x){all(!is.na(x))})
    sim <- r_truncated_multivariate_normal(n=1,
                                           alpha=c(as.matrix(alpha[ii,]))[iii],
                                           beta=c(as.matrix(beta[ii,]))[iii],
                                           gamma=t(tmp[iii,,drop=FALSE]),
                                           Gamma=x$Gamma,
                                           which="SZ",
                                           eps=1e-8,
                                           maxit=maxit)
    tmp <- rep(NA,sum(ii)*length(x$items)); tmp[iii] <- sim$S; residual[ii,] <- tmp
    Z[Z[[x$subject]]==i,2:(1+length(x$random.eff))] <- sim$Z
  }

  return(list(residuals=residual,ranef=Z))
}
