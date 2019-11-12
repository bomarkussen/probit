#' Simulate from a truncated multivariate normal distribution
#'
#' @description
#' Simulate \code{S in R^M} and \code{Z in R^d} given \code{S in (alpha,beta]$}, where \code{S = X + gamma^T Z} with \code{X ~ Normal(0,I_N)} and \code{Z ~ Normal(0,Gamma^T Gamma)}.
#' @param n Number of samples. Default: \code{n=1}.
#' @param alpha Vector of lower bounds for S component.
#' @param beta Vector of upper bounds for S component.
#' @param gamma Matrix of linear transformations on Z component.
#' @param Gamma Cholesky factor for variance on Z component.
#' @param which Text string (\code{"S"/"Z"/"SZ"}) deciding which componts are return. Default: \code{which="Z"}.
#' @param eps Convergence criteria for coupling from the past. Default: \code{eps=1e-12}.
#' @param maxit Maximal number of steps in Gibbs sampler. Default: \code{maxit=100}.
#' @return Simulations either as a matrix (if \code{which="S"} or \code{"Z"}), or as a list of matrices (if \code{which="SZ"}).
#' @examples
#' r_truncated_multivariate_normal(n=10,alpha=rep(-1,5),beta=rep(1,5),gamma=matrix(rnorm(10),2,5),Gamma=matrix(c(1,1,0,1),2,2))
#' @export
r_truncated_multivariate_normal <- function(n=1,alpha,beta,gamma,Gamma,which="Z",eps=1e-12,maxit=100) {
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
  if (eps <= 0) stop("eps must be strictly positive")
  if (maxit <= 0) stop("maxit must be strictly positive")

  # initialize factors
  d <- nrow(gamma)
  M <- ncol(gamma)
  Gamma.gamma <- Gamma%*%gamma
  delta0      <- solve(diag(nrow=d)+Gamma.gamma%*%t(Gamma.gamma),Gamma.gamma)
  delta       <- t(Gamma.gamma)%*%delta0
  Rt          <- t(chol(solve(diag(nrow=d)+Gamma.gamma%*%t(Gamma.gamma)))%*%Gamma)

  # Remark: delta is symmetric, cf. hack below

  # make matrices to contain the simulations
  if (is.element(which,c("S","SZ"))) S <- matrix(0,n,M)
  if (is.element(which,c("Z","SZ"))) Z <- matrix(0,n,d)

  # make simulations
  for (m in 1:n) {
    # initialize
    gammaT.RtV <- as.vector(t(gamma)%*%Rt%*%rnorm(d))
    f1 <- delta*alpha
    f2 <- delta*beta
    A  <- colSums(pmin(f1,f2)) + gammaT.RtV
    B  <- colSums(pmax(f1,f2)) + gammaT.RtV

    max.BA <- 0

    # Gibbs sampler
    for (iter in 1:maxit) {
      # random input
      U   <- runif(M)
      gammaT.RtV <- as.vector(t(gamma)%*%Rt%*%rnorm(d))

      # update of S
      x  <- pmin(A,B)
      GA <- x + qnorm(log_weighted_mean(U,pnorm(beta-x,log.p=TRUE),pnorm(alpha-x,log.p=TRUE)),log.p=TRUE)
      x  <- pmax(A,B)
      GB <- x + qnorm(log_weighted_mean(U,pnorm(beta-x,log.p=TRUE),pnorm(alpha-x,log.p=TRUE)),log.p=TRUE)

      # update of Z
      # remark: uses hack to do elementwise multiplication
      f1 <- delta*GA
      f2 <- delta*GB
      A  <- colSums(pmin(f1,f2)) + gammaT.RtV
      B  <- colSums(pmax(f1,f2)) + gammaT.RtV

      # convergence reached?
      if (max(abs(B-A))<eps) break
    }

    # keep track of coupling
    max.BA <- max(max.BA,max(abs(B-A)))

    # one additional step: find S
    U  <- runif(M)
    x  <- (A+B)/2
    Sm <- x + qnorm(log_weighted_mean(U,pnorm(beta-x,log.p=TRUE),pnorm(alpha-x,log.p=TRUE)),log.p=TRUE)

    # one additional step: find Z
    Zm <- as.vector(t(Gamma)%*%delta0%*%Sm + Rt%*%rnorm(d))

    # save results
    if (is.element(which,c("S","SZ"))) S[m,] <- Sm
    if (is.element(which,c("Z","SZ"))) Z[m,] <- Zm
  }

  # Did any iterations stop before coupling from the past was reached?
  if (max.BA>=eps) warning(paste("Coupling from the past not reached: max(|B-A|)=",signif(max.BA),". Consider increasing maxit.",sep=""))

  # return result
  if (which=="S") return(S)
  if (which=="Z") return(Z)
  if (which=="SZ") return(list(S=S,Z=Z))
}
