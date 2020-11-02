# !diagnostics off

#' Probit analysis with random effects
#'
#' @description
#' Separate analysis over items.
#' @param formula Model formula, where multivariate responses may be given additively on the left hand side. Responses must be ordered factors.
#' @param random List of formulas for the random effects. Models are fitted on first appearance observations.
#' @param subject Categorical variable encoding subjects.
#' @param dependence Text string (\code{"marginal" or "joint"}) deciding whether random effects are assumed independent or with a common joint normal distribution. Default: \code{dependence="marginal"}.
#' @param Gamma Choleskey factor for initial variance of random effects. If \code{Gamma=NULL} then initialized at identity matrix. Default: \code{Gamma=NULL}.
#' @param item.name Character string with name of generated variable identifying the items. Defaults to \code{item.name="item"}.
#' @param response.name Character string with name of generated variable containing the responses. Defaults to \code{response.name=NULL}, which corresponds to adding \code{".value"} to \code{item.name}.
#' @param data Date frame with data on the wide format.
#' @param data.long Possible additional long format data frame with item level explanatory variables. Presently not implemented!
#' @param B Number of simulations in minimization step. Default: \code{B=300}.
#' @param BB Number of simulations per subject in maximization step. Default: \code{BB=50}.
#' @param maxit Maximal number of minimization-maximization steps. Default: \code{steps=20}.
#' @param sig.level Significance level at which the iterative stochastic optimizations will be stopped. Defaults to \code{sig.level=0.60}.
#' @param verbose Numeric controlling amount of convergence diagnostics. Default: \code{verbose=0} corresponding to no output.
#' @note
#' A data frame must be provided, i.e. the \code{data} option is not optional. Variables not appearing in the \code{data}, but appearing in the \code{formula}, will be assumed to be random.
#' @return Object of class \code{probit}.
#' @details
#' If \code{Gamma=NULL} then predicted random effects are initialized at zero. Otherwise they are simulated from normal distribution with variance \code{Gamma^T Gamma}.
#' @export
probit <- function(formula,random,subject="id",dependence="marginal",Gamma=NULL,item.name="item",response.name=NULL,data,data.long=NULL,B=300,BB=50,maxit=20,sig.level=0.60,verbose=0) {
  # grab parameters ----

  # item and response name
  if ((!is.character(item.name))|(length(item.name)==0)) stop("item.name must be a non-vanishing character")
  if (is.null(response.name)) response.name <- paste(item.name,"value",sep=".")
  if ((!is.character(response.name))|(length(response.name)==0)) stop("response.name must be a non-vanishing character or NULL")

  # find items and make sanity check
  items <- all.vars(formula[[2]])
  if (!all(sapply(data[,items],function(x){is.numeric(x)|is.ordered(x)}))) stop("Response variables must be either numeric or ordered factor")
  items.interval <- items[sapply(data[,items],is.numeric)]
  items.ordinal  <- items[sapply(data[,items],is.ordered)]

  # find subjects
  subjects <- unique(data[[subject]])

  # names of random effects
  random.eff <- unlist(lapply(random,function(x){all.vars(x[[2]])}))
  q <- length(random.eff)
  if (q==0) stop("Present implementation assumes at least one random effect")

  # work with tibbles
  data <- as_tibble(data)

  # if random effects appear in data then there are removed
  if (length(intersect(names(data),random.eff))>0) {
    warning("Variables",paste(intersect(names(data),random.eff),collapse=", "),"are removed from data as they are used as random effects")
    data <- select(data,-any_of(random.eff))
  }

  # if item name appears in data then it is removed
  if (is.element(item.name,names(data))) {
    warning("Variable",item.name,"is removed from data it is used as item.name")
    data <- select(data,-any_of(item.name))
  }

  # if response name appears in data then it is removed
  if (is.element(response.name,names(data))) {
    warning("Variable",response.name,"is removed from data it is used as response.name")
    data <- select(data,-any_of(response.name))
  }

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

  # initial maximization step ----

  # mydata with random effects fixated at zero
  U <- as_tibble(cbind(subjects,matrix(0,length(subjects),q)),
                 .name_repair = "minimal")
  names(U) <- c(subject,random.eff)
  mydata <- full_join(as_tibble(data),U,by=subject)

  # simulate continuous responses for ordinal items
  for (i in items.ordinal) {
    # fit ordinal regression
    m.clm <- ordinal::clm(mydata[[i]]~1,link="probit",na.action = na.exclude)
    # simulated underlying normal and insert in data
    mydata[[i]] <- qnorm(c(0,pnorm(m.clm$alpha))[as.numeric(mydata[[i]])] + (diff(c(0,pnorm(m.clm$alpha),1))[as.numeric(mydata[[i]])])*runif(nrow(mydata)))
  }

  # linear regression
  mydata <- pivot_longer(mydata,all_of(items),names_to = item.name,values_to = response.name)
  m0 <- lm(eval(substitute(update(formula,y~.),list(y=as.name(response.name)))),
           data=mydata)

  # estimate sigma's
  mydata[!is.na(mydata[[response.name]]),response.name] <- residuals(m0)
  mydata <- pivot_wider(mydata,names_from = all_of(item.name), values_from = all_of(response.name))
  sigma2 <- vector("list",length(items.interval))
  names(sigma2) <- items.interval
  for (i in items.interval) {
    sigma2[[i]] <- mean((mydata[[i]]-mean(mydata[[i]],na.rm=TRUE))^2,na.rm=TRUE)
  }

  # predict with model and estimate threshold parameters
  # Remark: Reuses U from above
  eta        <- vector("list",length(items.ordinal))
  names(eta) <- items.ordinal
  mydata <- full_join(as_tibble(data),U,by=subject)
  for (i in items.ordinal) {
    tmp <- tibble(factor(rep(i,nrow(mydata)),levels=items)); names(tmp) <- item.name
    my.offset <- predict(m0,bind_cols(mydata,tmp))
    # estimate thresholds from ordinal regression
    m.clm <- ordinal::clm(mydata[[i]]~offset(my.offset),link="probit",na.action = na.exclude)
    eta[[i]] <- m.clm$alpha
  }

  # initialize mu and psi
  mu  <- matrix(0,length(subjects),q)
  psi <- matrix(0,length(subjects),q*(q+1)/2)
  for (s in 1:length(subjects)) psi[s,] <- Gamma[upper.tri(Gamma,diag=TRUE)]

  # random effects models
  lm.random <- vector("list",q); names(lm.random) <- random.eff
  data.short <- as_tibble(data) %>% group_by(!!as.name(subject)) %>% slice_head(n=1)
  U <- as.data.frame(cbind(subjects,mu)); names(U) <- c(subject,random.eff)
  data.short <- full_join(data.short,U,by=subject)
  for (i in 1:q) lm.random[[i]] <- lm(random[[i]],data=data.short)

  # means of random effects
  g <- matrix(0,length(subjects),q)
  for (i in 1:q) g[,i] <- predict(lm.random[[i]])


  # objective function for minimization step ----

  # linear parametrization matrix Q such that
  # 1) diagonal elements in parameter indexes = cumsum(1:q)
  # 2) Q Psi = matrix(Q%*%Psi,q,q)
  # 3) Q.mu  = matrix(tildeQ%*%mu,q,q*(q+1)/2)
  Q <- matrix(0,q*q,q*(q+1)/2)
  Q[matrix(1:(q*q),q,q)[upper.tri(matrix(0,q,q),diag=TRUE)],] <- diag(nrow=q*(q+1)/2)
  tildeQ <- matrix(aperm(array(Q,dim=c(q,q,q*(q+1)/2)),c(1,3,2)),q*q*(q+1)/2,q)

  # minimization function
  estimate.mu.psi <- function(s, ...) {
    # objective function
    F1 <- function(par,gradient=TRUE,hessian=TRUE) {
      # psi dimension
      r <- q*(q+1)/2

      # hack: replace zero diagonal by small number
      par[q+cumsum(1:q)] <- par[q+cumsum(1:q)] + (par[q+cumsum(1:q)]==0)*1e-4

      # unfold parameter
      mu.s  <- par[1:q]
      Psi.s <- matrix(Q%*%par[q+(1:r)],q,q)

      invPsi   <- solve(Psi.s)
      invPsi.Q <- invPsi%*%matrix(tildeQ,q,r*q)

      one.psi <- rep(0,r); one.psi[cumsum(1:q)] <- 1/par[q+cumsum(1:q)]

      # update data frame
      if (q > 1) {
        data.s[,random.eff] <-
          matrix(rep(matrix(mu.s,B,q,byrow=TRUE) + t(invPsi%*%U),
                     each=nrow(data.s)/B),nrow(data.s),q)
      } else {
        data.s[,random.eff] <-
          rep(matrix(mu.s,B,q,byrow=TRUE) + t(invPsi%*%U),
              each=nrow(data.s)/B)
      }
      data.s.long <- pivot_longer(data.s,all_of(items), names_to = item.name, values_to = response.name)

      # compute sum over observations
      inner.sum <- colSums(matrix(
        mapply(function(name,value,f){
          ifelse(is.element(name,items.interval),
                 log(sigma2[[name]])/2+((value-f)^2)/(2*sigma2[[name]]),
                 -log(pnorm((eta[[name]])[value]-f)-
                        pnorm(c(-Inf,eta[[name]])[value]-f)))},
          data.s.long[[item.name]],data.s.long[[response.name]],
          predict(m0,newdata=data.s.long)),
        nrow(data.s.long)/B,B),na.rm=TRUE)

      # hack: set infinite inner sums to a large number
      inner.sum[is.infinite(inner.sum)] <- 1e4

      # compute F1
      res <- sum(log(abs(par[q+cumsum(1:q)]))) - q/2 - log(det(Gamma)) +
        sum((Gamma%*%(mu.s-g[s,]))^2)/2 + sum(c(Gamma%*%invPsi)^2)/2 +
        sum(inner.sum)/B

      # compute gradient?
      if (gradient) {
        Gamma.invPsi.Q.invPsi <-
          matrix(aperm(array(matrix(Gamma%*%invPsi.Q,q*r,q)%*%invPsi,dim=c(q,r,q)),
                       c(1,3,2)),q*q,r)
        attr(res,"gradient") <-
          c(t(Gamma)%*%Gamma%*%(mu.s-g[s,]),
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
    U <- matrix(rnorm(q*B),q,B); rownames(U) <- random.eff

    # set-up data matrix
    data.s <- merge(filter(data,(!!as.name(subject))==subjects[s]),t(U),by=NULL)
    for (i in items.ordinal) data.s[[i]] <- as.numeric(data.s[[i]])

    # minimize F1
    res <- nlm(F1,c(mu[s,],psi[s,]),check.analyticals = FALSE)

    # return
    return(c(res$minimum,res$estimate))
  }

  # MM-loop ----
  code <- 1
  for (iter in 1:maxit) {
    # minimization step ----

    # minimizations allowing for parallization via future::plan()
    res <- future_map(1:length(subjects),estimate.mu.psi,.progress=(verbose>1),.options = furrr_options(seed = TRUE))
    if (verbose > 1) cat("\n")

    F1.last <- sapply(res,function(x){x[1]})
    mu  <- t(sapply(res,function(x){x[1+(1:q)]}))
    psi <- t(sapply(res,function(x){x[-(1:(q+1))]}))

    # convergence diagnostics
    if (iter==1) {
      if (verbose>0) cat("iteration",iter,": F1=",sum(F1.last),"\n")
    } else {
      pval <- t.test(F1.last-F1.prev,alternative="less")$p.value
      if (verbose>0) cat("iteration",iter,": F1=",sum(F1.last),", change=",sum(F1.last-F1.prev),", p-value(no decrease)=",pval,"\n")
    }

    # maximization step ----

    # set-up data matrix with random input
    U <- as_tibble(cbind(rep(subjects,each=BB),matrix(0,length(subjects)*BB,q)),
                   .name_repair = "minimal")
    names(U) <- c(subject,random.eff)
    for (s in 1:length(subjects)) {
      U[U[[subject]]==subjects[s],-1] <- t(mu[s,] + solve(matrix(Q%*%psi[s,],q,q),matrix(rnorm(BB*q),q,BB)))
    }
    mydata <- full_join(as_tibble(data),U,by=subject)

    # predict with previous model and simulate responses for ordinal variables
    for (i in items.ordinal) {
      tmp <- tibble(factor(rep(i,nrow(mydata)),levels=items)); names(tmp) <- item.name
      my.offset <- predict(m0,bind_cols(mydata,tmp))
      # fit ordinal regression
      m.clm <- ordinal::clm(mydata[[i]]~offset(my.offset),
                            link="probit",na.action = na.exclude)
      # simulated underlying normal and insert in mydata
      mydata[[i]] <- my.offset + qnorm(
        pnorm(c(-Inf,m.clm$alpha)[as.numeric(mydata[[i]])]-my.offset) +
          (pnorm(c(m.clm$alpha,Inf)[as.numeric(mydata[[i]])]-my.offset) - pnorm(c(-Inf,m.clm$alpha)[as.numeric(mydata[[i]])]-my.offset))*runif(nrow(mydata))
        )
    }

    # linear regression
    mydata <- pivot_longer(mydata,all_of(items),names_to = item.name, values_to = response.name)
    m0 <- lm(eval(substitute(update(formula,y~.),list(y=as.name(response.name)))),
             data=mydata)

    # estimate sigma's
    mydata[!is.na(mydata[[response.name]]),response.name] <- residuals(m0)
    mydata <- pivot_wider(mydata,names_from = all_of(item.name), values_from = all_of(response.name))
    for (i in items.interval) {
      sigma2[[i]] <- mean((mydata[[i]]-mean(mydata[[i]],na.rm=TRUE))^2,na.rm=TRUE)
    }

    # predict with previous model and update threshold parameters
    # Remark: Reuses random input U from above
    mydata <- full_join(as_tibble(data),U,by=subject)
    for (i in items.ordinal) {
      tmp <- tibble(factor(rep(i,nrow(mydata)),levels=items)); names(tmp) <- item.name
      my.offset <- predict(m0,bind_cols(mydata,tmp))
      # estimate thresholds from ordinal regression
      m.clm <- ordinal::clm(mydata[[i]]~offset(my.offset),
                            link="probit",na.action = na.exclude)
      eta[[i]] <- m.clm$alpha
    }

    # random effects models
    data.short <- as_tibble(data) %>% group_by(!!as.name(subject)) %>% slice_head(n=1)
    U <- as.data.frame(cbind(subjects,mu)); names(U) <- c(subject,random.eff)
    data.short <- full_join(data.short,U,by=subject)
    for (i in 1:q) lm.random[[i]] <- lm(random[[i]],data=data.short)

    # means of random effects
    for (i in 1:q) g[,i] <- predict(lm.random[[i]])

    # Estimate Gamma
    hat.var <- matrix(rowMeans(
      apply(mu-g,1,function(x){x%*%t(x)}) +
      apply(psi,1,function(x){solve(t(matrix(Q%*%x,q,q))%*%matrix(Q%*%x,q,q))})
      ),q,q)
    if (dependence!="joint") {
      hat.var[upper.tri(hat.var)] <- 0
      hat.var[lower.tri(hat.var)] <- 0
    }
    Gamma <- chol(solve(hat.var))
    rownames(Gamma) <- colnames(Gamma) <- random.eff

    # break MM-loop?
    if ((iter>1) && (pval>sig.level)) {
      code <- 0
      break
    }

    # update and end EM-loop
    F1.prev <- F1.last
  }

  # return probit-object ----
  return(structure(list(item.name=item.name,items.interval=items.interval,items.ordinal=items.ordinal,
                        ranef=random.eff,dependence=dependence,
                        m.fixed=m0,sigma2=sigma2,eta=eta,m.random=lm.random,Gamma=Gamma,
                        mu=mu,psi=psi,
                        B=B,BB=BB,F1=F1.last,pvalue=pval,iter=iter,code=code,
                        data=data),class="probit"))
}
