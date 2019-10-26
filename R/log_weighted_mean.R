#' Weighted mean on logarithmic scale
#'
#' @param U Vector of weights in \code{[0,1]}.
#' @param log.a Vector of logarithms of the first term.
#' @param log.b Vector of logarithms of the second term.
#' @return Logarithm of \code{U*a+(1-U)*b}.
#' @examples
#' log_weighted_mean(U=c(1,0.5,0),log.a=c(-10,12,42),log.b=c(60,2,-10))
#' log_weighted_mean(U=0.4,log.a=-2.88e-21,log.b=-1.38e-16)
#' log_weighted_mean(U=0.295,log.a=-5.748586e-23,log.b=-1.199389e-18)
#' log_weighted_mean(U=0.5751,log.a=-6.414215e-22,log.b=-2.806906e-17)
#' @export
log_weighted_mean <- function(U,log.a,log.b) {
  tmp <- ifelse((U>0)&(log(U/(1-U))>log.b-log.a),
                log1p((1-U)/U*exp(log.b-log.a))+log1p(-1+U)+log.a,
                log1p(U/(1-U)*exp(log.a-log.b))+log1p(-U)+log.b)
  pmax(pmin(log.a,log.b),pmin(pmax(log.a,log.b),tmp))
}
