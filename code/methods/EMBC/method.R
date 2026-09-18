library(parallel)
library(glmnet)
set.seed(123)
sed = sample(1:1e8,1000)
get_true_weibull <- function(a0, a1, b) {
  if (a0 <= 0 || a1 <= 0 || b <= 0) {
    stop("a0, a1, and b must all be positive.")
  }
  
  k   <- b
  beta  <- b * (a0^(-b) - a1^(-b))
  beta0 <- b * log(a0 / a1) + (a0^(-b) - a1^(-b))
  AUC <- a1^k / (a0^k + a1^k)
  
  opt_cut <- as.numeric(inv_boxcox(-beta0 / beta, k))
  
  se = 1 - pweibull(opt_cut, shape = b, scale = a1)
  sp = pweibull(opt_cut, shape = b, scale = a0)
  J = se + sp - 1
  
  out = c(beta0=beta0,beta=beta,k=k,AUC=AUC,opt_cut=opt_cut,J=J,se=se,sp=sp)
  out
}
gendata_weibull <- function(n0=100, n1=100, pi0=0.9,pi1=0.9,a0 = 1/2, a1 = 2, b = 1.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (a0 <= 0 || a1 <= 0 || b <= 0) {
    stop("s0, s1, and k must all be positive.")
  }
  
  n00=rbinom(1,n0,pi0)
  n11=rbinom(1,n1,pi1)
  X = c(rweibull(n00, shape = b, scale = a0),rweibull(n0-n00, shape = b, scale = a1))
  Y = c(rweibull(n11, shape = b, scale = a1),rweibull(n1-n11, shape = b, scale = a0))
  
  true_param <- get_true_weibull(a0 = a0, a1 = a1, b = b)
  
  list(T = c(X,Y), R = c(rep(0, n0),rep(1,n1)), G = c(rep(0,n00),rep(1,n0-n00),rep(1,n11),rep(0,n1-n11)),
       true_param = true_param
  )
}


gendata_lognormal <- function(n0 = 100, n1 = 100, pi0 = 0.9, pi1 = 0.9, a0 = 0, a1 = 1, b = 1,seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (b <= 0) {
    stop("b must be positive.")
  }

  n00 = rbinom(1, n0, pi0)
  n11 = rbinom(1, n1, pi1)
  X = c(
    rlnorm(n00, meanlog = a0, sdlog = sqrt(b)),
    rlnorm(n0 - n00, meanlog = a1, sdlog = sqrt(b))
  )
  Y = c(
    rlnorm(n11, meanlog = a1, sdlog = sqrt(b)),
    rlnorm(n1 - n11, meanlog = a0, sdlog = sqrt(b))
  )

  true_param <- get_true_lognormal(a0 = a0, a1 = a1, b = b)

  list(T = c(X,Y), R = c(rep(0, n0),rep(1,n1)), G = c(rep(0,n00),rep(1,n0-n00),rep(1,n11),rep(0,n1-n11)),
    true_param = true_param
  )
}

gendata_beta <- function(n0 = 100, n1 = 100, pi0 = 0.9, pi1 = 0.9,a0 = 1.5, a1 = 3, b = 3,seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (a0 <= 0 || b <= 0 || a1 <= 0 ) {
    stop("a0, b0, a1, and b1 must all be positive.")
  }

  n00 = rbinom(1, n0, pi0)
  n11 = rbinom(1, n1, pi1)
  X = c(
    rbeta(n00, shape1 = a0, shape2 = b),
    rbeta(n0 - n00, shape1 = a1, shape2 = b)
  )
  Y = c(
    rbeta(n11, shape1 = a1, shape2 = b),
    rbeta(n1 - n11, shape1 = a0, shape2 = b)
  )

  true_param <- get_true_beta(a0 = a0, a1 = a1, b = b)

  list(T = c(X,Y), R = c(rep(0, n0),rep(1,n1)), G = c(rep(0,n00),rep(1,n0-n00),rep(1,n11),rep(0,n1-n11)),
       true_param = true_param
  )
}
gendata_gamma <- function(n0 = 100, n1 = 100, pi0 = 0.9, pi1 = 0.9, a0 = 1.5, a1 = 3, b = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (a0 <= 0 || b <= 0 || a1 <= 0) {
    stop("a0, b, and a1 must all be positive.")
  }

  n00 = rbinom(1, n0, pi0)
  n11 = rbinom(1, n1, pi1)
  X = c(
    rgamma(n00, shape = a0, rate = b),
    rgamma(n0 - n00, shape = a1, rate = b)
  )
  Y = c(
    rgamma(n11, shape = a1, rate = b),
    rgamma(n1 - n11, shape = a0, rate = b)
  )

  true_param <- get_true_gamma(a0 = a0, a1 = a1, b = b)

  list(T = c(X,Y), R = c(rep(0, n0),rep(1,n1)), G = c(rep(0,n00),rep(1,n0-n00),rep(1,n11),rep(0,n1-n11)),
       true_param = true_param
  )
}

gendata_gamma2 <- function(n0 = 100, n1 = 100, pi0 = 0.9, pi1 = 0.9, a = 1, b0 = 0.5, b1 = 1.73, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (a <= 0 || b0 <= 0 || b1 <= 0) {
    stop("a, b0, and b1 must all be positive.")
  }
  if (b0 == b1) {
    stop("b0 and b1 must be different (unequal-rate gamma). Use gendata_gamma for equal-rate case.")
  }

  n00 = rbinom(1, n0, pi0)
  n11 = rbinom(1, n1, pi1)
  X = c(
    rgamma(n00, shape = a, rate = b0),
    rgamma(n0 - n00, shape = a, rate = b1)
  )
  Y = c(
    rgamma(n11, shape = a, rate = b1),
    rgamma(n1 - n11, shape = a, rate = b0)
  )

  true_param <- get_true_gamma2(a = a, b0 = b0, b1 = b1)

  list(T = c(X,Y), R = c(rep(0, n0),rep(1,n1)), G = c(rep(0,n00),rep(1,n0-n00),rep(1,n11),rep(0,n1-n11)),
       true_param = true_param
  )
}
get_true_lognormal = function(a0, a1, b){
  beta0 = (a0^2-a1^2)/2/b
  beta = (a1-a0)/b
  k=0
  AUC = pnorm((a1-a0)/sqrt(2*b))
  opt_cut = as.numeric(exp(-beta0/beta))
  J = 2*pnorm(abs(a0-a1)/2/sqrt(b))-1
  se = 1-plnorm(opt_cut,a1,sqrt(b)) 
  sp = plnorm(opt_cut, a0, sqrt(b))
  c(beta0=beta0,beta=beta,k=k,AUC=AUC,opt_cut=opt_cut,J=J,se=se,sp=sp)
}
get_true_gamma = function(a0, a1, b){
  beta0 = (a1-a0)*log(b)-log(gamma(a1))+log(gamma(a0))
  beta = (a1-a0)
  k=0
  AUC = pbeta(0.5, shape1 = a0, shape2 = a1)
  opt_cut = as.numeric(exp(-beta0/beta))
  J = abs(pgamma(opt_cut, shape = a0,rate = b)-pgamma(opt_cut, shape = a1,rate = b))
  se = 1-pgamma(opt_cut, shape = a1,rate = b)
  sp = pgamma(opt_cut, shape = a0,rate = b)
  c(beta0=beta0,beta=beta,k=k,AUC=AUC,opt_cut=opt_cut,J=J,se=se,sp=sp)
}
get_true_beta = function(a0, a1, b){
  beta0 = log(gamma(a1+b))-log(gamma(a1))+log(gamma(a0)) - log(gamma(a0+b))
  beta = (a1-a0)
  k=0
  M = 1e5
  x0 = rbeta(M, shape1 = a0, shape2 = b)
  x1 = rbeta(M, shape1 = a1, shape2 = b)
  AUC =   mean(x1 > x0)
  opt_cut = as.numeric(exp(-beta0/beta))
  J = abs(pbeta(opt_cut, shape1 = a0,shape2 = b)-pbeta(opt_cut, shape1 = a1,shape2 = b))
  se = 1-pbeta(opt_cut, shape1 = a1,shape2 = b)
  sp = pbeta(opt_cut, shape1 = a0,shape2 = b)
  c(beta0=beta0,beta=beta,k=k,AUC=AUC,opt_cut=opt_cut,J=J,se=se,sp=sp)
}

get_true_gamma2 = function(a, b0, b1){
  beta = b0-b1
  beta0 = a*log(b1/b0)+beta
  k=1
  AUC = pbeta(b1 / (b0 + b1), shape1 = a, shape2 = a)
  opt_cut = as.numeric(inv_boxcox(-beta0/beta, k))
  if (b1 > b0) {
    se = pgamma(opt_cut, shape = a, rate = b1)
    sp = 1-pgamma(opt_cut, shape = a, rate = b0)
  } else {
    se = 1-pgamma(opt_cut, shape = a, rate = b1)
    sp = pgamma(opt_cut, shape = a, rate = b0)
  }
  J = se + sp - 1
  c(beta0=beta0,beta=beta,k=k,AUC=AUC,opt_cut=opt_cut,J=J,se=se,sp=sp)
}

log1pexp <- function(x) {
  out <- numeric(length(x))
  hi <- is.finite(x) & x > 0
  out[hi] <- x[hi] + log1p(exp(-x[hi]))
  out[!hi] <- log1p(exp(x[!hi]))
  out
}

boxcox = function(x, lam){
  if (abs(lam) < 1e-8) return(log(x))
  (x^lam-1)/lam
}
inv_boxcox = function(x,lam){
  if (abs(lam) < 1e-8) return(exp(x))
  (x*lam + 1)^(1/lam)
}



get_rbmse = function(info,true){
  k = length(true)
  rb = sapply(1:k, FUN = function(i) {
    if(true[i] ==0) mean((info[i,] - true[i]))
    mean((info[i,] - true[i])/true[i])
    })
  mse = sapply(1:k, FUN = function(i) mean((info[i,] - true[i])^2))
  rbind(true,rb,mse)
}


get_est = function(fit,data, EM = FALSE){
  temp = fit
  if (EM) {
    temp = fit$best
    pi0 = temp$pi0
    pi1 = temp$pi1
  }
  T = data$T
  R = data$R
  n0 = sum(R == 0)
  n1 = sum(R == 1)
  n = n0+n1
  alpha = temp$alpha
  beta = temp$beta
  k = temp$k
  pi = temp$pi
  qi = pi*exp(alpha+beta*boxcox(T, k))
  pi=pi/sum(pi);qi=qi/sum(qi)
  F0hat = function(x,alpha,beta,T,pi) sum(pi*I((alpha+beta*boxcox(T,k))<x))
  F1hat = function(x,alpha,beta,T,qi) sum(qi*I((alpha+beta*boxcox(T,k))<x))
  
  auc1 = sum(sapply(1:n,FUN = function(i) qi[i]*F0hat((alpha+beta*boxcox(T[i],k)),alpha,beta,T,pi)))
  auc2 = 1-sum(sapply(1:n,FUN = function(i) pi[i]*F1hat((alpha+beta*boxcox(T[i],k)),alpha,beta,T,qi)))
  auc = (auc1+auc2)/2
  
  opt_cut = inv_boxcox(-alpha/beta,k)
  
  J = abs(F0hat((alpha+beta*boxcox(opt_cut,k)),alpha,beta,T,pi)-F1hat((alpha+beta*boxcox(opt_cut,k)),alpha,beta,T,qi))
  
  se = 1-F1hat((alpha+beta*boxcox(opt_cut,k)),alpha,beta,T,qi)
  sp = F0hat((alpha+beta*boxcox(opt_cut,k)),alpha,beta,T,pi)
  
  if (EM){
    return(c(alpha = alpha,beta = beta, k=k, AUC = auc, opt_cut = opt_cut, J=J,se = se, sp=sp, pi0=pi0, pi1=pi1))
  } else{
    return(  c(alpha = alpha,beta = beta, k=k, AUC = auc, opt_cut = opt_cut, J=J,se = se, sp=sp))
  }
}


pl = function(data,w0, w1, k, beta_pen = 1/2, pi_pen = 1/2){
  t = data$T
  R = data$R
  n0 = sum(R == 0)
  n1 = sum(R == 1)
  n=n0+n1
  x = t[data$R == 0] 
  y = t[data$R == 1]
  lam_beta = n^-beta_pen 
  lam_pi0 = n0^-pi_pen
  lam_pi1 = n1^-pi_pen
  
  pi0 = (sum(w0)+lam_pi0)/(n0+lam_pi0)
  pi1 = (sum(w1)+lam_pi1)/(n1+lam_pi1)
  lambda = (sum(1-w0) + sum(w1))/n
  
  res = c(rep(1, n), rep(0, n))
  cov = t(rbind(rep(0, 2*n), c(boxcox(x, k), boxcox(y, k), boxcox(x, k), boxcox(y, k))))
  # Keep both weighted classes strictly represented for glmnet. During EM,
  # responsibilities can underflow to zero, which makes glmnet reject the
  # weighted response with error 9001 (null probability too close to 1).
  weights <- pmax(c(1-w0, w1, w0, 1-w1), 1e-4)
  
  fit_temp = glmnet(y = res, x = cov, family = "binomial", weights = weights, lambda = 2*lam_beta/n, standardize = FALSE)
  alphastar = coef(fit_temp)[1]
  beta = coef(fit_temp)[3]
  alpha = alphastar - log(n*lambda/(n-n*lambda))
  pi = (n+(sum(1-w0)+sum(w1))*(exp(alpha+beta*boxcox(t, k))-1))^-1
  
  pi_part = sum(log(pi))
  x_part = sum(sapply(x, FUN = function(x) log(pi0+(1-pi0)*exp(alpha+beta*boxcox(x, k)))))
  y_part = sum(sapply(y, FUN = function(y) log(1-pi1+pi1*exp(alpha+beta*boxcox(y, k)))))
  loglik = pi_part + x_part + y_part
  ploglik = pi_part + x_part + y_part + lam_pi0*log(pi0)+lam_pi1*log(pi1)-lam_beta*beta%*%beta
  list(alpha = alpha, beta = beta, pi0=pi0, pi1 = pi1, loglik = loglik, ploglik = ploglik,pi=pi)
}



EMBC_single = function(data, beta_pen = 1/2, pi_pen = 1/2, pi0, pi1, alpha, beta, k, k_interval = c(-2, 2), 
                       tol = 1e-6, maxite = 1000){
  T = data$T
  R = data$R
  n0 = sum(R == 0)
  n1 = sum(R == 1)
  n=n0+n1
  x = T[data$R == 0] 
  y = T[data$R == 1]
  lam_beta = n^-beta_pen 
  lam_pi0 = n0^-pi_pen
  lam_pi1 = n1^-pi_pen
  loglikold = 100; logliknew = 0
  err = loglikold-logliknew 
  ite = 0
  
  while (err > tol & ite < maxite){
    eta0 = alpha + beta*boxcox(x, k)
    eta1 = alpha + beta*boxcox(y, k)
    w0 = if (pi0 == 0) rep(0, length(eta0)) else if (pi0 == 1) rep(1, length(eta0)) else plogis(qlogis(pi0) - eta0)
    w1 = if (pi1 == 0) rep(0, length(eta1)) else if (pi1 == 1) rep(1, length(eta1)) else plogis(qlogis(pi1) + eta1)
    
    k = optimize(
      f = function(k) -pl(data, w0, w1, k,beta_pen = beta_pen, pi_pen = pi_pen)$ploglik,
      interval = k_interval
    )$minimum
    
    out = pl(data, w0, w1, k,beta_pen = beta_pen, pi_pen = pi_pen)
    alpha = out$alpha
    beta = out$beta
    pi0 = out$pi0
    pi1 = out$pi1 
    logliknew = out$ploglik
    loglik = out$loglik
    pi = out$pi
    
    ite = ite + 1
    err = abs(logliknew - loglikold)
    loglikold = logliknew 
  }
  list(pi0=pi0, pi1 = pi1,alpha = alpha, beta = beta, k = k,loglik = loglik,ploglik = logliknew, pi = pi, numite = ite)
}

EMBC = function(data, beta_pen = 1/2, pi_pen = 1/2, k_interval = c(-2, 2), 
                tol = 1e-6, maxite = 1000){
  T = data$T
  R = data$R
  pi0 = seq(0.6, 0.9, by = 0.1)
  pi1 = seq(0.6, 0.9, by = 0.1)
  
  fit = glm(R~T, family=binomial)
  alpha = fit$coefficients[1]
  beta = fit$coefficients[2]
  
  fit = list();loglik = c()
  for (i in 1:4){
    pi0_temp = pi0[i];pi1_temp = pi1[i]
    k = runif(10, -1.5, 1.5)
    for(j in 1:5){
      k_temp = k[j]
      alpha_temp = alpha + runif(1, -3, 3)
      beta_temp = beta + runif(1, -3, 3)
      fit_temp = tryCatch(
        EMBC_single(data, beta_pen = beta_pen, pi_pen = pi_pen,
                    pi0 = pi0_temp, pi1 = pi1_temp,
                    alpha = alpha_temp, beta = beta_temp, k = k_temp,
                    k_interval = k_interval, tol = tol, maxite = maxite),
        error = function(e) NULL
      )
      if (is.null(fit_temp) || length(fit_temp$loglik) != 1L ||
          !is.finite(fit_temp$loglik)) next
      fit = c(fit, list(fit_temp))
      loglik = c(loglik, fit_temp$loglik)
    }
  }

  if (length(fit) == 0L) stop("All EMBC random starts failed.")
  
  ind= which.max(loglik)
  best = fit[[ind]]  
  list(best = best, fit= fit,loglik = loglik)
}





pl_fixpi = function(data,pi0, pi1,w0, w1, k, beta_pen = 1/2){
  t = data$T
  R = data$R
  n0 = sum(R == 0)
  n1 = sum(R == 1)
  n=n0+n1
  x = t[data$R == 0] 
  y = t[data$R == 1]
  lam_beta = n^-beta_pen 
  
  lambda = (sum(1-w0) + sum(w1))/n
  
  res = c(rep(1, n), rep(0, n))
  cov = t(rbind(rep(0, 2*n), c(boxcox(x, k), boxcox(y, k), boxcox(x, k), boxcox(y, k))))
  # Keep both weighted classes strictly represented for glmnet. During EM,
  # responsibilities can underflow to zero, which makes glmnet reject the
  # weighted response with error 9001 (null probability too close to 1).
  weights <- pmax(c(1-w0, w1, w0, 1-w1), 1e-4)
  
  fit_temp = glmnet(y = res, x = cov, family = "binomial", weights = weights, lambda = 2*lam_beta/n, standardize = FALSE)
  alphastar = coef(fit_temp)[1]
  beta = coef(fit_temp)[3]
  alpha = alphastar - log(n*lambda/(n-n*lambda))
  pi = (n+(sum(1-w0)+sum(w1))*(exp(alpha+beta*boxcox(t, k))-1))^-1
  
  pi_part = sum(log(pi))
  x_part = sum(sapply(x, FUN = function(x) log(pi0+(1-pi0)*exp(alpha+beta*boxcox(x, k)))))
  y_part = sum(sapply(y, FUN = function(y) log(1-pi1+pi1*exp(alpha+beta*boxcox(y, k)))))
  loglik = pi_part + x_part + y_part
  ploglik = pi_part + x_part + y_part - lam_beta*beta%*%beta
  list(alpha = alpha, beta = beta, pi0=pi0, pi1 = pi1, loglik = loglik, ploglik = ploglik,pi=pi)
}



EMBC_single_fixpi = function(data, beta_pen = 1/2, pi0, pi1, alpha, beta, k, k_interval = c(-2, 2), 
                       tol = 1e-6, maxite = 1000){
  T = data$T
  R = data$R
  n0 = sum(R == 0)
  n1 = sum(R == 1)
  n=n0+n1
  x = T[data$R == 0] 
  y = T[data$R == 1]
  lam_beta = n^-beta_pen 
  loglikold = 100; logliknew = 0
  err = loglikold-logliknew 
  ite = 0
  
  while (err > tol & ite < maxite){
    eta0 = alpha + beta*boxcox(x, k)
    eta1 = alpha + beta*boxcox(y, k)
    w0 = if (pi0 == 0) rep(0, length(eta0)) else if (pi0 == 1) rep(1, length(eta0)) else plogis(qlogis(pi0) - eta0)
    w1 = if (pi1 == 0) rep(0, length(eta1)) else if (pi1 == 1) rep(1, length(eta1)) else plogis(qlogis(pi1) + eta1)
    
    k = optimize(
      f = function(k) -pl_fixpi(data, pi0, pi1, w0, w1, k,beta_pen = beta_pen)$ploglik,
      interval = k_interval
    )$minimum
    
    out = pl_fixpi(data,pi0,pi1, w0, w1, k,beta_pen = beta_pen)
    alpha = out$alpha
    beta = out$beta
    logliknew = out$ploglik
    loglik = out$loglik
    pi = out$pi
    
    ite = ite + 1
    err = abs(logliknew - loglikold)
    loglikold = logliknew 
  }
  list(pi0=pi0, pi1 = pi1,alpha = alpha, beta = beta, k = k,loglik = loglik,ploglik = logliknew, pi = pi, numite = ite)
}

EMBC_fixpi = function(data, pi0, pi1, beta_pen = 1/2, k_interval = c(-2, 2), 
                tol = 1e-6, maxite = 1000){
  T = data$T
  R = data$R
  
  fit = glm(R~T, family=binomial)
  alpha = fit$coefficients[1]
  beta = fit$coefficients[2]
  
  fit = list();loglik = c()
  k = runif(20, -1.5, 1.5)
  for(j in 1:20){
      k_temp = k[j]
      alpha_temp = alpha + runif(1, -3, 3)
      beta_temp = beta + runif(1, -3, 3)
      fit_temp = tryCatch(
        EMBC_single_fixpi(data, beta_pen = beta_pen, pi0 = pi0, pi1 = pi1,
                          alpha = alpha_temp, beta = beta_temp, k = k_temp,
                          k_interval = k_interval, tol = tol, maxite = maxite),
        error = function(e) NULL
      )
      if (is.null(fit_temp) || length(fit_temp$loglik) != 1L ||
          !is.finite(fit_temp$loglik)) next
      fit = c(fit, list(fit_temp))
      loglik = c(loglik, fit_temp$loglik)
    }

  if (length(fit) == 0L) stop("All EMBC random starts failed.")
  
  
  ind= which.max(loglik)
  best = fit[[ind]]  
  list(best = best, fit= fit,loglik = loglik)
}


# Fit the EMBC model with known class proportions and a fixed Box-Cox
# parameter. The unified bootstrap uses this after estimating k once from the
# original sample.
EMBC_single_fixpi_fixk = function(data, beta_pen = 1/2, pi0, pi1,
                                  alpha, beta, k, tol = 1e-6,
                                  maxite = 1000) {
  T = data$T
  R = data$R
  x = T[R == 0]
  y = T[R == 1]
  loglikold = 100
  logliknew = 0
  err = abs(loglikold - logliknew)
  ite = 0

  while (err > tol & ite < maxite) {
    eta0 = alpha + beta*boxcox(x, k)
    eta1 = alpha + beta*boxcox(y, k)
    w0 = if (pi0 == 0) rep(0, length(eta0)) else if (pi0 == 1) rep(1, length(eta0)) else plogis(qlogis(pi0) - eta0)
    w1 = if (pi1 == 0) rep(0, length(eta1)) else if (pi1 == 1) rep(1, length(eta1)) else plogis(qlogis(pi1) + eta1)

    out = pl_fixpi(data, pi0, pi1, w0, w1, k,
                   beta_pen = beta_pen)
    alpha = out$alpha
    beta = out$beta
    logliknew = out$ploglik
    loglik = out$loglik
    pi = out$pi

    ite = ite + 1
    err = abs(logliknew - loglikold)
    loglikold = logliknew
  }

  list(pi0 = pi0, pi1 = pi1, alpha = alpha, beta = beta, k = k,
       loglik = loglik, ploglik = logliknew, pi = pi, numite = ite)
}


EMBC_fixpi_fixk = function(data, pi0, pi1, k, beta_pen = 1/2,
                           tol = 1e-6, maxite = 1000) {
  if (length(k) != 1 || !is.finite(k)) {
    stop("k must be one finite numeric value.")
  }

  T = data$T
  R = data$R
  initial_fit = glm(R ~ T, family = binomial)
  alpha = initial_fit$coefficients[1]
  beta = initial_fit$coefficients[2]

  fit = list()
  loglik = c()
  for (j in 1:20) {
    fit_temp = tryCatch(
      EMBC_single_fixpi_fixk(
        data, beta_pen = beta_pen, pi0 = pi0, pi1 = pi1,
        alpha = alpha + runif(1, -3, 3),
        beta = beta + runif(1, -3, 3), k = k,
        tol = tol, maxite = maxite
      ),
      error = function(e) NULL
    )
    if (is.null(fit_temp) || length(fit_temp$loglik) != 1L ||
        !is.finite(fit_temp$loglik)) next
    fit = c(fit, list(fit_temp))
    loglik = c(loglik, fit_temp$loglik)
  }

  if (length(fit) == 0L) stop("All EMBC random starts failed.")

  ind = which.max(loglik)
  best = fit[[ind]]
  list(best = best, fit = fit, loglik = loglik)
}
