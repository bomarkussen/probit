#' Gibbs sampler for a truncated multivariate normal distribution
#'
#' @description
#' Simulate \code{S in R^M} and \code{Z in R^d} given \code{S in (alpha,beta]$}, where \code{S = X + gamma^T Z} with \code{X ~ Normal(0,I_N)} and \code{Z ~ Normal(0,Gamma^T Gamma)}.
#' @param n Number of samples. Default: \code{n=1}.
#' @param alpha Vector of lower bounds for S component.
#' @param beta Vector of upper bounds for S component.
#' @param gamma Matrix of linear transformations on Z component.
#' @param Gamma Cholesky factor for variance on Z component.
#' @param which Text string (\code{"S"/"Z"/"SZ"}) deciding which componts are return. Default: \code{which="Z"}.
#' @param steps Number of steps in Gibbs sampler. Default: \code{steps=20}.
#' @return Simulations either as a matrix (if \code{which="S"} or \code{"Z"}), or as a list of matrices (if \code{which="SZ"}).
#' @examples
#' Gibbs_truncated_multivariate_normal(n=10,alpha=rep(-1,5),beta=rep(1,5),gamma=matrix(rnorm(10),2,5),Gamma=matrix(c(1,1,0,1),2,2))
#' @export
Gibbs_truncated_multivariate_normal <- function(n=1,alpha,beta,gamma,Gamma,which="Z",steps=20) {
  #
  # with n=10,000 the above example took 26.72 seconds
  #
  # sanity checks
  # 1. alpha and beta must be vectors of the same length, M
  # 2. alpha must be coordinatewise strictly less then beta
  # 3. gamma must be a matrix with ncols = M
  # 4. Gamma must square matrix with same nrows as gamma, d
  if (!is.vector(alpha)) stop("alpha must be a vector")
  if (!is.vector(beta)) stop("beta must be a vector")
  if (!is.matrix(gamma)) stop("gamma must be a matrix")
  if (!is.matrix(Gamma)) stop("Gamma must be a matrix")
  if (length(alpha)!=length(beta)) stop("alpha and beta must be of the same length")
  if (length(alpha)!=ncol(gamma)) stop("gamma must have number of columns equal to length of alpha")
  if (nrow(Gamma)!=ncol(Gamma)) stop("Gamma must be a square matrix")
  if (nrow(Gamma)!=nrow(gamma)) stop("gamma and Gamma must have same number of rows")
  if (any(alpha>=beta)) stop("alpha must be coordinate-wise strictly less than beta")
  if (!is.element(which,c("S","Z","SZ"))) stop("which must be either S, Z or SZ")
  if (steps <= 0) stop("steps must be strictly positive")

  # initialize factors
  d <- nrow(gamma)
  M <- ncol(gamma)
  Gamma.gamma <- Gamma%*%gamma
  delta0      <- solve(diag(nrow=d)+Gamma.gamma%*%t(Gamma.gamma),Gamma.gamma)
  delta       <- t(Gamma.gamma)%*%delta0
  Rt          <- t(chol(solve(diag(nrow=d)+Gamma.gamma%*%t(Gamma.gamma)))%*%Gamma)

  # Remark: delta is symmetric, cf. hack below

  # make simulations
  S <- replicate(n,{
    # initialize
    G <- (alpha+beta)/2

    # Gibbs sampler
    for (iter in 1:steps) {
      # random input
      U <- runif(M)
      gammaT.RtV <- as.vector(t(gamma)%*%Rt%*%rnorm(d))

      # update of Z
      Z  <- colSums(delta*G) + gammaT.RtV

      # update of S
      G <- Z + qnorm(log_weighted_mean(U,pnorm(beta-Z,log.p=TRUE),pnorm(alpha-Z,log.p=TRUE)),log.p=TRUE)
    }

    # return result
    G
  })

  # one additional step: find Z
  Z <- apply(matrix(S,M,n),2,function(x){as.vector(t(Gamma)%*%delta0%*%x + Rt%*%rnorm(d))})

  # return result
  if (which=="S") return(S=t(S))
  if (which=="Z") return(Z=t(Z))
  if (which=="SZ") return(list(S=t(S),Z=t(Z)))
}
