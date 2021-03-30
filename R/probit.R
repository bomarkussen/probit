# !diagnostics off

#' Probit analysis with random effects
#'
#' @description
#' Separate analysis over items.
#' @param fixed Model formula for the fixed effect, where multivariate responses may be given additively on the left hand side. Responses must be ordered factors.
#' @param random List of formulas for the random effects. Models are fitted on first appearance observations.
#' @param subject.name Character string with name of categorical variable encoding subjects.
#' @param dependence Text string (\code{"marginal" or "joint"}) deciding whether random effects are assumed independent or with a common joint normal distribution. Default: \code{dependence="marginal"}.
#' @param Gamma Choleskey factor for initial variance of random effects. If \code{Gamma=NULL} then initialized at identity matrix. Default: \code{Gamma=NULL}.
#' @param item.name Character string with name of generated variable identifying the items. Defaults to \code{item.name="item"}.
#' @param response.name Character string with name of generated variable containing the responses. Defaults to \code{response.name=NULL}, which corresponds to adding \code{".value"} to \code{item.name}.
#' @param data Date frame with data on the wide format.
#' @param data.long Possible additional long format data frame with item level explanatory variables. Presently not implemented!
#' @param mu Matrix (no subjects, q) of initial estimates for mu.
#' @param psi Matrix (no subjects, q*(q+1)/2) of initial estimates for psi.
#' @param B Number of simulations in minimization step. Default: \code{B=300}.
#' @param BB Number of simulations per subject in maximization step. Default: \code{BB=50}.
#' @param maxit Maximal number of minimization-maximization steps. Default: \code{maxit=20}.
#' @param sig.level Significance level at which the iterative stochastic optimizations will be stopped. Defaults to \code{sig.level=0.60}.
#' @param verbose Numeric controlling amount of convergence diagnostics. Default: \code{verbose=0}. Level \code{1} gives progress bar and convergence statistics. Level \code{2} also gives diagnostics for maximization step. Level \code{3} also gives graphical diagnostics for minimization step.
#' @note
#' A data frame must be provided, i.e. the \code{data} option is not optional. Variables that
#' coincide with random effects, item identifier or internal name for item response will be removed
#' from \code{data}, and an warning will be issued.
#'
#' The minimization step is implemented via \code{\link[furrr]{future_map}}. This implies that the user may
#' activate parallel computations by calling \code{\link[future]{plan}}, e.g. \code{future::plan("multisession", workers = 4)}.
#'
#' @return \code{\link{probit-class}} object.
#' @details
#'
#' @export
probit <- function(fixed,random,subject.name="id",dependence="marginal",Gamma=NULL,item.name="item",response.name=NULL,data,data.long=NULL,mu=NULL,psi=NULL,B=300,BB=50,maxit=20,sig.level=0.60,verbose=0) {
  # sanity check and grab parameters ----

  # data.long has not yet been implemented
  if (!is.null(data.long)) warning("Additional long format data frame has not yet been implemented")

  # item and response name
  if ((!is.character(item.name))|(length(item.name)==0)) stop("item.name must be a non-vanishing character")
  if (is.null(response.name)) response.name <- paste(item.name,"value",sep=".")
  if ((!is.character(response.name))|(length(response.name)==0)) stop("response.name must be a non-vanishing character or NULL")

  # TO DO: weight name
  weight.name <- paste(item.name,"weight",sep=".")

  # find items and make sanity check
  items <- all.vars(fixed[[2]])
  if (!all(sapply(data[,items],function(x){is.numeric(x)|is.factor(x)}))) stop("Response variables must be either numeric or factor")
  items.interval <- items[sapply(data[,items],is.numeric)]
  items.ordinal  <- items[sapply(data[,items],is.factor)]

  # names of random effects
  random.eff <- sapply(random,function(x){all.vars(x[[2]])})
  q <- length(random.eff)
  if (q==0) stop("Present implementation assumes at least one random effect")

  # work with tibbles
  data <- as_tibble(data)

  # take levels of factors and recode them as numeric
  ordinal.levels <- as.list(rep(NA,length(items.ordinal)))
  names(ordinal.levels) <- items.ordinal
  for (i in items.ordinal) {
    ordinal.levels[[i]] <- levels(data[[i]])
    data[[i]] <- as.numeric(data[[i]])
  }

  # if random effects appear in data then these are removed
  if (length(intersect(names(data),random.eff))>0) {
    warning("Variables ",paste(intersect(names(data),random.eff),collapse=", ")," are removed from data as they are used as random effects")
    data <- select(data,-any_of(random.eff))
  }

  # if item.name appears in data then it is removed
  if (is.element(item.name,names(data))) {
    warning("Variable ",item.name," is removed from data as it is used as item.name")
    data <- select(data,-any_of(item.name))
  }

  # if response.name appears in data then it is removed
  if (is.element(response.name,names(data))) {
    warning("Variable ",response.name," is removed from data as it is used as response.name")
    data <- select(data,-any_of(response.name))
  }

  # if weight.name appears in data then it is removed
  if (is.element(weight.name,names(data))) {
    warning("Variable ",weight.name," is removed from data as it is used as weight.name")
    data <- select(data,-any_of(weight.name))
  }

  # insert column to contain weight.name
  tmp <- tibble(weight.name=1); names(tmp) <- weight.name
  data <- full_join(data,tmp,by=character())

  # remove observations with missing explanatory variables
  data[[weight.name]] <- as.numeric(complete.cases(select(data,any_of(setdiff(
    c(all.vars(fixed),unlist(sapply(random,function(x){all.vars(x)}))),
    c(items,random.eff))))))
  if (sum(1-data[[weight.name]])>0) {
    warning("Remove ",sum(1-data[[weight.name]])," observations with non-complete explanatory variables. NOTE: This may interact badly with update().")
    data <- filter(data,!!as.name(weight.name)==1)
  }

  # find subjects
  subjects <- unique(data[[subject.name]])

  # fix dependence
  if (dependence!="marginal") dependence <- "joint"

  # Cholesky factor for variance of random effects
  if (is.null(Gamma)) {
    # if not prespecified, then initialize as the identity matrix
    Gamma <- diag(nrow=q)
  } else {
    # check format of Gamma
    if (nrow(Gamma)!=ncol(Gamma)) stop("Gamma must be NULL or a square matrix")
    if (nrow(Gamma)!=q) stop("Gamma must have number of columns equal to the number of random effects")
  }

  # Mean values in Gaussian approximation of conditional distribution of random effects
  if (is.null(mu)) {
    mu  <- matrix(0,length(subjects),q)
  } else {
    if (!is.matrix(mu)) stop("If specified, then mu must be a matrix")
    if (nrow(mu)!=length(subjects)) stop("Number of rows in mu must match number of subjects")
    if (ncol(mu)!=q) stop("Number of columns in mu must match number of random effects")
  }

  # Cholesky factor for precision in Gaussian approximation of conditional distribution of random effects
  if (is.null(psi)) {
    psi <- matrix(0,length(subjects),q*(q+1)/2)
    for (s in 1:length(subjects)) psi[s,] <- Gamma[upper.tri(Gamma,diag=TRUE)]
  } else {
    if (!is.matrix(psi)) stop("If specified, then psi must be a matrix")
    if (nrow(psi)!=length(subjects)) stop("Number of rows in psi must match number of subjects")
    if (ncol(psi)!=(q*(q+1)/2)) stop("Number of columns in psi must be q*(q+1)/2, where q is number of random effects")
  }

  # return result from Minimization-Maximization iterations
  return(MM_probit(maxit,sig.level,verbose,
                   fixed,response.name,weight.name,
                   item.name,items.interval,items.ordinal,ordinal.levels,
                   subject.name,random,dependence,
                   m.fixed=NULL,sigma2=NULL,eta=NULL,m.random=NULL,Gamma=Gamma,
                   mu,psi,
                   B,BB,
                   data,
                   estimate.models = TRUE))
}


# help function for probit() and update.probit() ----

# log-function: log(1-exp(-a))
# Martin Mächer, log1mexp-note, ETH Zurich, October 2012.
# log1mexp <- function(a) {ifelse(a<0.693,log(-expm1(-a)),log1p(-exp(-a)))}

# log normal probability
log_pnorm <- function(a,b) {
  mysign <- -sign(pmin(a,b))
  log.b  <- pnorm(pmax(mysign*a,mysign*b),log.p=TRUE)
  log.a  <- pnorm(pmin(mysign*a,mysign*b),log.p=TRUE)
  log.b + ifelse(log.b-log.a<0.693,log(-expm1(log.a-log.b)),log1p(-exp(log.a-log.b)))
}

# probit analysis criterion function
my.clm <- function(par,y,ff,logF) {
  #  tmp <- cumsum(c(0,rep(par[1],length(ordinal.levels[[i]])-2)))
  tmp <- cumsum(c(par[2],abs(par[-(1:2)])+0.001))
  -sum(logF((c(tmp,Inf)[y]-ff)/sqrt(abs(par[1])+0.001),
            (c(-Inf,tmp)[y]-ff)/sqrt(abs(par[1])+0.001)))
}



# Maximization-Minimization for probit model

MM_probit <- function(maxit,sig.level,verbose,
                      fixed,response.name,weight.name,
                      item.name,items.interval,items.ordinal,ordinal.levels,
                      subject.name,random,dependence,
                      m.fixed,sigma2,eta,m.random,Gamma,
                      mu,psi,
                      B,BB,
                      data,
                      estimate.models=TRUE) {
  # if estimate.models=FALSE, then only one iteration with a minimization step
  # and without maximization step is done.
  if (!estimate.models) {
    maxit <- 1
    logL  <- 0
  }

  # grab parameters
  random.eff <- unlist(lapply(random,function(x){all.vars(x[[2]])}))
  q          <- length(random.eff)
  subjects   <- unique(data[[subject.name]])
  items      <- c(items.interval,items.ordinal)

  # if NULL, then initialize sigma2
  if (is.null(sigma2)) {
    sigma2 <- as.list(rep(1,length(items)))
    names(sigma2) <- items
  }

  # if NULL, then initialize eta
  if (is.null(eta)) {
    eta        <- vector("list",length(items))
    names(eta) <- items
    for (i in items.ordinal) {
      tmp <- pmax(1,table(factor(data[[i]],levels=ordinal.levels[[i]])))
      eta[[i]] <- qnorm(cumsum(tmp[-length(tmp)])/sum(tmp))
    }
  }

  # if NULL, the initialize m.fixed
  if (is.null(m.fixed)) {
    # mydata with predicted random effects
    U <- as_tibble(cbind(subjects,mu), .name_repair = "minimal")
    names(U) <- c(subject.name,random.eff)
    mydata <- full_join(data,U,by=subject.name)

    # simulate continuous responses for ordinal items
    for (i in items.ordinal) {
      # simulated underlying normal and insert in data
      a <- (c(-Inf,eta[[i]])[mydata[[i]]] - 0)/sqrt(sigma2[[i]])
      b <- (c(eta[[i]],Inf)[mydata[[i]]] - 0)/sqrt(sigma2[[i]])
      flip <- -sign(a)
      a <- pnorm(flip*a)
      b <- pnorm(flip*b)
      mydata[[i]] <- 0 + sqrt(sigma2[[i]])*flip*qnorm(pmin(a,b) + abs(b-a)*runif(nrow(mydata)))
    }

    # linear regression
    mydata  <- pivot_longer(mydata,all_of(items),names_to = item.name,values_to = response.name)
    mydata[[item.name]] <- factor(mydata[[item.name]],levels=items)
#    mydata <- mydata[!is.na(mydata[[response.name]]),]
    m.fixed <- biglm::biglm(eval(substitute(update(fixed,y~.),list(y=as.name(response.name)))),
                            data=mydata)
#    m.fixed <- lm(eval(substitute(update(fixed,y~.),list(y=as.name(response.name)))),
#                            data=mydata)
  }

  # if NULL, then initialize m.random
  if (is.null(m.random)) {
    m.random <- vector("list",length(random.eff))
    names(m.random) <- random.eff
  }

  # linear parametrization matrix Q such that
  # 1) diagonal elements in parameter indexes = cumsum(1:q)
  # 2) Q Psi = matrix(Q%*%Psi,q,q)
  # 3) Q.mu  = matrix(tildeQ%*%mu,q,q*(q+1)/2)
  Q <- matrix(0,q*q,q*(q+1)/2)
  Q[matrix(1:(q*q),q,q)[upper.tri(matrix(0,q,q),diag=TRUE)],] <- diag(nrow=q*(q+1)/2)
  tildeQ <- matrix(aperm(array(Q,dim=c(q,q,q*(q+1)/2)),c(1,3,2)),q*q*(q+1)/2,q)

  # maximization function for ordinal responses ----
  estimate.eta <- function(i, ...) {
    # manual encoding of likelihood function
    y <- mydata[[i]][!is.na(mydata[[i]])]
    tmp <- tibble(factor(i,levels=items))
    names(tmp) <- item.name
    my.offset  <- predict_slim(m.fixed.new,full_join(mydata,tmp,by=character()))[!is.na(mydata[[i]])]
    my.res <- optim(par=c(sigma2[[i]]-0.001,eta[[i]][1],diff(eta[[i]])-0.001),
                    fn=my.clm,gr=NULL,
                    y=y,ff=my.offset,logF=log_pnorm)
    # thresholds
    my.levels     <- ordinal.levels[[i]]
    my.eta        <- cumsum(c(my.res$par[2],abs(my.res$par[-(1:2)])+0.001))
    names(my.eta) <- paste(my.levels[-length(my.levels)],my.levels[-1],sep="|")
    # variance
    my.sigma2     <- abs(my.res$par[1])+0.001

    # return result
    return(list(eta=my.eta,sigma2=my.sigma2,logL= -my.res$value/BB, entropy=my.res$value/(length(y)*BB)))
  }

  # minimization function ----
  estimate.mu.psi <- function(s, ...) {
    # objective function
    F1 <- function(par,gradient=!TRUE,hessian=!TRUE) {
      # psi dimension
      r <- q*(q+1)/2

      # hack: replace zero diagonal by small number
      par[q+cumsum(1:q)] <- par[q+cumsum(1:q)] + (par[q+cumsum(1:q)]==0)*1e-8

      # unfold parameter
      mu.s  <- par[1:q]
      Psi.s <- matrix(Q%*%par[q+(1:r)],q,q)

      invPsi   <- solve(Psi.s)
      invPsi.Q <- invPsi%*%matrix(tildeQ,q,r*q)

      one.psi <- rep(0,r); one.psi[cumsum(1:q)] <- 1/par[q+cumsum(1:q)]

      # update data frame
      data.s[,random.eff] <- matrix(mu.s,nrow(data.s),q,byrow=TRUE) + U.s%*%t(invPsi)
      data.s.long <- pivot_longer(data.s,all_of(items), names_to = item.name, values_to = response.name)
      data.s.long <- data.s.long[!is.na(data.s.long[[response.name]]),]
      # first use item.name as a factor
      data.s.long[[item.name]] <- factor(data.s.long[[item.name]],levels=items)
      tmp <- predict_slim(m.fixed,newdata=data.s.long)
      # thereafter recode as numeric
      data.s.long[[item.name]] <- as.numeric(data.s.long[[item.name]])

      # compute sum over observations
      inner.sum <- colSums(matrix(
        mapply(function(name,value,f){
          ifelse(is.element(items[name],items.interval),
                 log(sigma2[[name]])/2+((value-f)^2)/(2*sigma2[[name]]),
                 -log_pnorm((c(eta[[name]],Inf)[value]-f)/sqrt(sigma2[[name]]),
                            (c(-Inf,eta[[name]])[value]-f)/sqrt(sigma2[[name]]))
                 )},
          data.s.long[[item.name]],data.s.long[[response.name]],tmp),
        nrow(data.s.long)/B,B),na.rm=TRUE)

      # hack: set infinite inner sums to a large number
      if (any(is.infinite(inner.sum))) {
        warning(paste0(sum(is.infinite(inner.sum)),"infinites in inner sum"))
        inner.sum[is.infinite(inner.sum)] <- 1e6
      }

      # compute F1
      res <- sum(log(abs(par[q+cumsum(1:q)]))) - q/2 - log(det(Gamma)) +
        sum((Gamma%*%(mu.s-mean.Z[s,]))^2)/2 + sum(c(Gamma%*%invPsi)^2)/2 +
        sum(inner.sum)/B

      # compute gradient?
      if (gradient) {
        Gamma.invPsi.Q.invPsi <-
          matrix(aperm(array(matrix(Gamma%*%invPsi.Q,q*r,q)%*%invPsi,dim=c(q,r,q)),
                       c(1,3,2)),q*q,r)
        attr(res,"gradient") <-
          c(t(Gamma)%*%Gamma%*%(mu.s-mean.Z[s,]),
            one.psi - t(Gamma.invPsi.Q.invPsi)%*%c(Gamma%*%invPsi)) +
          rbind(t(Psi.s)%*%U,
                apply(U,2,function(x){one.psi - t(matrix(tildeQ%*%invPsi%*%x,q,r))%*%x}))%*%inner.sum/B
      }

      # compute hessian?
      if (gradient & hessian) {
        # analytic part
        hes <- matrix(0,q+r,q+r)
        hes[1:q,1:q] <- t(Gamma)%*%Gamma
        tmp <- array(0,dim=c(r,r))
        for (i in 1:r) for (j in 1:r) tmp[i,j] <-
          sum(c(Gamma%*%invPsi)*c(Gamma%*%invPsi%*%
                                    (matrix(Q[,i],q,q)%*%invPsi%*%matrix(Q[,j],q,q)+
                                       matrix(Q[,j],q,q)%*%invPsi%*%matrix(Q[,i],q,q))%*%invPsi))
        diag(tmp)[cumsum(1:q)] <- diag(tmp)[cumsum(1:q)] - 1/(par[q+cumsum(1:q)]^2)
        hes[q+(1:r),q+(1:r)] <- tmp +
          t(Gamma.invPsi.Q.invPsi)%*%Gamma.invPsi.Q.invPsi
        # stochastic part
        tmp.mu.mu <- matrix(apply(U,2,function(x){x%*%t(x)-diag(nrow=q)})%*%inner.sum/B,q,q)
        tmp.psi.mu <- matrix(apply(U,2,function(x){
          one.psi%*%t(x)%*%Psi.s +
            t(matrix(tildeQ%*%invPsi%*%x,q,r))%*%(diag(nrow=q)-x%*%t(x))%*%Psi.s +
            matrix(matrix(t(Q),r*q,q)%*%x,r,q)})%*%inner.sum/B,r,q)
        tmp.psi.psi <- matrix(apply(U,2,function(x){
          t(matrix(tildeQ%*%invPsi%*%x,q,r))%*%(diag(nrow=q)-x%*%t(x))%*%matrix(tildeQ%*%invPsi%*%x,q,r) -
            one.psi%*%t(x)%*%matrix(tildeQ%*%invPsi%*%x,q,r) -
            t(matrix(tildeQ%*%invPsi%*%x,q,r))%*%x%*%t(one.psi)})%*%inner.sum/B,r,r)
        # combine result
        hes[1:q,1:q] <- hes[1:q,1:q] + tmp.mu.mu
        hes[q+(1:r),1:q] <- hes[q+(1:r),1:q] + tmp.psi.mu
        hes[1:q,q+(1:r)] <- hes[1:q,q+(1:r)] + t(tmp.psi.mu)
        hes[q+(1:r),q+(1:r)] <- hes[q+(1:r),q+(1:r)] + tmp.psi.psi
        attr(res,"hessian") <- hes
      }

      # return result
      return(res)
    }

    # sample random input
#    U <- as_tibble(cbind(rep(subjects[s],each=B),t(matrix(UU[,1:B,s],q,B))),
#                   .name_repair = "minimal")
    U <- as_tibble(cbind(rep(subjects[s],each=B),t(matrix(rnorm(q*B),q,B))),
                   .name_repair = "minimal")
    names(U) <- c(subject.name,random.eff)

    # set-up data matrix and random input for F1
    data.s <- full_join(filter(data,(!!as.name(subject.name))==subjects[s]),
                        U,by=subject.name)
    U.s <- as.matrix(data.s[,random.eff])
    U   <- t(as.matrix(U[,random.eff]))

    # minimize F1 if there are any observations
    if (any(!is.na(c(as.matrix(data.s[,c(items.interval,items.ordinal)]))))) {
#      res <- nlm(F1,c(mu[s,],psi[s,]),check.analyticals = FALSE)
      res <- optim(c(mu[s,],psi[s,]),F1)
      res <- list(minimum=res$value,estimate=res$par)
    } else {
      res <- list(minimum=0,estimate=c(mu[s,],psi[s,]))
    }

    # return
    return(c(res$minimum,res$estimate))
  }

  # Hack: Fix sample of normal distributions
  #UU <- array(rnorm(q*max(B,BB)*length(subjects)),dim=c(q,max(B,BB),length(subjects)))

  # MM-loop ----
  code <- 1
  pval <- NA
  for (iter in 1:sum(maxit)) {
    # maximization step ----

    # estimate model parameters?
    if (estimate.models) {
      # estimate random effects models
      U <- as_tibble(cbind(subjects,mu),.name_repair = "minimal")
      names(U) <- c(subject.name,random.eff)
      data.short <- data %>% group_by(!!as.name(subject.name)) %>% slice_head(n=1)
      data.short <- full_join(data.short,U,by=subject.name)
      for (i in 1:q) {
        m.random[[i]] <- lm(random[[i]],data=data.short)
        data.short[[random.eff[i]]] <- residuals(m.random[[i]])
      }

      # estimate Gamma
      #      hat.var <- matrix(rowMeans(matrix(
      #        apply(mu-mean.Z,1,function(x){x%*%t(x)}) +
      #          apply(psi,1,function(x){solve(t(matrix(Q%*%x,q,q))%*%matrix(Q%*%x,q,q))})
      #        ,q*q,nrow(mu))),q,q)
      hat.var <- var(data.short[,random.eff])*(1-(1/nrow(data.short))) +
        matrix(rowMeans(matrix(apply(psi,1,function(x){solve(t(matrix(Q%*%x,q,q))%*%matrix(Q%*%x,q,q))}),
                               q*q,nrow(psi))),q,q)
      if (dependence=="marginal") {
        hat.var[upper.tri(hat.var)] <- 0
        hat.var[lower.tri(hat.var)] <- 0
      }
      Gamma <- chol(solve(hat.var))
      rownames(Gamma) <- colnames(Gamma) <- random.eff

      # reinitiate mu via two-component maxit?
      if ((iter==maxit[1]) & (iter<sum(maxit))) {
        for (i in 1:q) mu[,i] <- predict(m.random[[i]])
      }

      # fixed effects estimation below

      # set-up data matrix with random input
      U <- as_tibble(cbind(rep(subjects,each=BB),matrix(0,length(subjects)*BB,q)),
                     .name_repair = "minimal")
      names(U) <- c(subject.name,random.eff)

      # set-up random input for ordinal responses
      my.runif <- vector("list",length(items.ordinal))
      names(my.runif) <- items.ordinal
      for (i in items.ordinal) my.runif[[i]] <- runif(nrow(data)*BB)

      # iterations of maximization step (at least 1 step done, at most 10 steps done)
      logL <- -Inf
      for (M.iter in 1:10) {
        # simulate random effects and set-up data frame
        for (s in 1:length(subjects)) {
          U[U[[subject.name]]==subjects[s],-1] <- t(mu[s,] + solve(matrix(Q%*%psi[s,],q,q),matrix(rnorm(q*BB),q,BB)))
          #U[U[[subject.name]]==subjects[s],-1] <- t(mu[s,] + solve(matrix(Q%*%psi[s,],q,q),matrix(UU[,1:BB,s],q,BB)))
        }
        mydata <- full_join(data,U,by=subject.name)

        # predict with previous model and simulate responses for ordinal variables
        for (i in items.ordinal) {
          tmp <- tibble(factor(i,levels=items))
          names(tmp) <- item.name
          my.offset  <- predict_slim(m.fixed,full_join(mydata,tmp,by=character()))

          # simulated underlying normal and insert in data
          a <- (c(-Inf,eta[[i]])[mydata[[i]]] - my.offset)/sqrt(sigma2[[i]])
          b <- (c(eta[[i]],Inf)[mydata[[i]]] - my.offset)/sqrt(sigma2[[i]])
          flip <- -sign(a)
          a <- pnorm(flip*a)
          b <- pnorm(flip*b)
          mydata[[i]] <- my.offset + sqrt(sigma2[[i]])*flip*qnorm(pmin(a,b) + abs(b-a)*my.runif[[i]])
        }

        # linear regression
        mydata <- pivot_longer(mydata,all_of(items),names_to = item.name, values_to = response.name)
        mydata[[item.name]]   <- factor(mydata[[item.name]],levels=items)
        mydata[[weight.name]] <- 1/unlist(sigma2)[as.numeric(mydata[[item.name]])]
        m.fixed.new <- biglm::biglm(eval(substitute(update(formula(fixed),y~.),list(y=as.name(response.name)))),
                                    weight = eval(substitute(formula(~w),list(w=as.name(weight.name)))),
                                    data = filter(mydata,!is.na(!!as.name(response.name))))
#        m.fixed <- lm(eval(substitute(update(formula(fixed),y~.),list(y=as.name(response.name)))),
#                                weight = weighted.OLS, data=mydata)

        # initiate log(likelihood)
        if (verbose>1) cat("Step",M.iter,":")
        logL.new   <- 0
        sigma2.new <- sigma2
        eta.new    <- eta

        # estimate sigma's for interval responses using 'mydata' from above
        mydata[[response.name]] <- mydata[[response.name]] - predict_slim(m.fixed.new,newdata=mydata)
        for (i in items.interval) {
          tmp <- (mydata[[response.name]])[mydata[[item.name]]==i]
          sigma2.new[[i]] <- mean(tmp^2,na.rm=TRUE)
          # add log(likelihood)
          logL.new <- logL.new - sum(!is.na(tmp))*(log(sigma2.new[[i]])+1)/(2*BB)
          if (verbose>1) cat(sqrt(sigma2.new[[i]]),",")
        }

        # estimate eta's, possibly in parallel session
#        future::plan("multisession", workers = 4, gc=TRUE)
        mydata <- full_join(data,U,by=subject.name)
        res    <- furrr::future_map(items.ordinal,estimate.eta,.options = furrr::furrr_options(seed = TRUE))
        names(res) <- items.ordinal
#        future:::ClusterRegistry("stop")

        # extract results for ordinal regressions
        logL.new <- logL.new + sum(sapply(res,function(x){x$logL}))
        eta.new[items.ordinal]    <- sapply(res,function(x){x$eta})
        sigma2.new[items.ordinal] <- sapply(res,function(x){x$sigma2})
        if (verbose>1) cat(sapply(res,function(x){x$entropy}))

        # end Miter
        if (verbose>1) cat(": logL=",logL.new,"\n")
        if (logL > logL.new) break
        logL    <- logL.new
        m.fixed <- m.fixed.new
        sigma2  <- sigma2.new
        eta     <- eta.new
      }

      # end model estimations
      rm(data.short,mydata,logL.new,m.fixed.new,sigma2.new,eta.new)
    }

    # minimization step ----

    # means of random effects
    mean.Z <- matrix(0,length(subjects),q)
    for (i in 1:q) {
      mean.Z[,i] <- predict_slim(m.random[[i]],
                                 slice_head(group_by(data,!!as.name(subject.name)),n=1))
    }

    # minimizations allowing for parallization via future::plan()
#    future::plan("multisession", workers = 4, gc=TRUE)
    res <- furrr::future_map(1:length(subjects),estimate.mu.psi,.progress=(verbose > 1),.options = furrr::furrr_options(seed = TRUE))
    F1.new <- sapply(res,function(x){x[1]})
#    future:::ClusterRegistry("stop")

    # convergence diagnostics
    if (verbose > 1) cat("\n")
    if (iter==1) {
      if (verbose > 0) cat("iteration",iter,": logL=",logL,", F1=",sum(F1.new),"\n")
    } else {
      pval <- t.test(F1.new-F1.best,alternative="less")$p.value
      if (verbose > 0) cat("iteration",iter,": logL=",logL,", F1=",sum(F1.new),", change=",sum(F1.new-F1.best),", p-value(no decrease)=",pval,"\n")
    }

    if (verbose>2) {
      if (iter==1) {
        print(ggplot(data.frame(x=F1.new),aes(x=x)) + geom_density(fill="lightgreen") +
                xlab(paste0("Iteration 1: sum(F1)=",round(sum(F1.new),2))))
      } else {
        print(ggExtra::ggMarginal(ggplot(data.frame(x=F1.best,y=F1.new),aes(x=x,y=y)) +
                                    geom_point() + geom_abline(intercept=0, slope=1, col="red") +
                                    xlab(paste0("Iteration ",iter-1,": sum(F1)=",round(sum(F1.best),2))) +
                                    ylab(paste0("Iteration ",iter,": sum(F1)=",round(sum(F1.new),2))) +
                                    xlim(min(F1.best,F1.new),max(F1.best,F1.new)) +
                                    ylim(min(F1.best,F1.new),max(F1.best,F1.new)),
                                  type="density",fill="lightgreen"))

      }
    }

    # break MM-loop?
    if ((iter>1) && (pval>sig.level)) {
      code <- 0
      break
    }

    # update parameters in minimization-step
    F1.best <- F1.new
    mu  <- t(matrix(sapply(res,function(x){x[1+(1:q)]}),q,length(res)))
    psi <- t(matrix(sapply(res,function(x){x[-(1:(q+1))]}),q*(q+1)/2,length(res)))

    # end MM-loop
  }

  # TO DO: In principle logL should be updated if MM-loop reaches maxiter

  # naming
  colnames(mu) <- random.eff

  # recode ordinal variables
  for (i in items.ordinal) {
    data[[i]] <- ordered((ordinal.levels[[i]])[data[[i]]],levels=ordinal.levels[[i]])
  }

  # return probit-object ----
  return(structure(list(fixed=fixed,response.name=response.name,weight.name=weight.name,
                        item.name=item.name,items.interval=items.interval,items.ordinal=items.ordinal,
                        ordinal.levels=ordinal.levels,
                        subject.name=subject.name,random=random,dependence=dependence,
                        m.fixed=m.fixed,sigma2=sigma2,eta=eta,m.random=m.random,Gamma=Gamma,
                        mu=mu,psi=psi,
                        B=B,BB=BB,logL=logL,F1=F1.best,pvalue=pval,iter=iter,code=code,
                        data=data),class="probit"))
}
