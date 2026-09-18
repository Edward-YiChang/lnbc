###EM Algorithm (given tuning parameter nu)
DRest0.EM <- function(dat,pi0,pi1,
                      nu=1,k=30,ord=4,ord.pen=2,
                      maxit=500,thres=1e-4,
                      knots=NULL,ini.G=NULL){
  ##################
  ###INPUTS:
  #dat: dataframe consists of R and biomarker 
  #pi0=P(G=0|R=0); pi1=P(G=1|R=1)
  #nu: smoothness penalty tuning parameter
  #k: number of B-spline basis
  #ord, ord.pen: the orders of B-splines and the penalty
  #maxit: maximum number of iterations
  #thres: threshold value for stopping
  #knots: user specified knots. 
  #       The number should be k+ord+1; the middle k-ord+1 must include all data points!!!
  #       "q": automatically generate knots according to
  #ini.G: two components specifying initial pseudo-responses for Groups R=0,1
  #       or G sequence with length matches sample size

  ###OUTPUTS: a list of:
  #fin.all: final gamObject output
  #obj.df: information for each iteration: J value, improvement of J, convergence status
  #conv: convergence status for EM
  #knots: knots used
  #vshift: log(pi1*m+(1-pi0)*n)-log((1-pi1)*m+pi0*n), useful for final h estimation
  #EL.tunes: tuning parameter for EL approach
  #weight.Dis, weight.Dis.est: weights vectors for empirical distribution
  #obsLL: observed log-likelihood of the final iteration
  ##################
  cov <- dat$biomarker
  id0 <- which(dat$R==0);id1 <- which(dat$R==1)
  n <- length(id0);m <- length(id1)
  vshift <- log(pi1*m+(1-pi0)*n)-log((1-pi1)*m+pi0*n)
  lam.prop <- (pi1*m+(1-pi0)*n)/(n+m)
  conv <- vector(mode = "logical", length =0)
  Jval <-  Jascend <- numeric(length = 0)
  dh.abs.mean <- dh.abs.max <- dh.abs.re <- numeric(length = 0)
  
  ###initialization:
  if(is.null(ini.G)){
    G.work <- ifelse(dat$R==1,pi1,1-pi0) #initial value for unknown G_i
  }else if(length(ini.G)==2){
    G.work <- ifelse(dat$R==1,ini.G[2],ini.G[1])
  }else{
    G.work <- ini.G
  }
  G.l <- list();G.l[[1]] <- G.work
  
  if(!is.null(knots)){
    if(knots=="q"){knots <- gen.q.knots(cov,k=k,ord=ord)}
  }
  
  gamout0 <- gam(G.work~s(cov,bs='bs',k=k,m=c(ord,ord.pen)),
                 family = binomial,sp=nu,
                 knots = list(cov = knots),
                 control = list(scalePenalty=F))#initial h estimation
  conv[1] <- gamout0$converged
  knots.his <- list(gamout0$smooth[[1]]$knots)
  #CAUSION: I use Ghat as response directly; the weights are constant
  #compute initial J value (multiplied by n+m)
  Jcomputation <- obs.LL.pen(id0,id1,gamout0,nu,pi0,pi1,lam.prop,vshift)
  Jval[1] <- Jcomputation$J
  exph.work <- Jcomputation$exp.signal
  #h evaluations:
  eval.interval <- quantile(dat$biomarker,probs=c(0.05,0.95))
  eval.pts <- seq(from=eval.interval[1],to=eval.interval[2],length.out=30)
  ht.est <- h.eval(gamout0,eval.pts,vshift)
  
  ###EM algorithm
  for(r in 1:maxit){
    #E-step
    G.work[id0] <- (1-pi0)*exph.work[id0]/((1-pi0)*exph.work[id0]+pi0)
    G.work[id1] <- pi1*exph.work[id1]/(pi1*exph.work[id1]+1-pi1)
    G.work <- pmax(G.work,0)
    G.work <- pmin(G.work,1)
    G.l[[r+1]] <- G.work
    #M-step
    gamout <- gam(G.work~s(cov,bs='bs',k=k,m=c(ord,ord.pen)),
                  family = binomial,sp=nu,
                  knots = list(cov = knots),
                  control = list(scalePenalty=F))
    conv[r+1] <- gamout$converged
    knots.his <- unique(c(knots.his,list(gamout$smooth[[1]]$knots)))
    
    Jcomputation <- obs.LL.pen(id0,id1,gamout,nu,pi0,pi1,lam.prop,vshift)
    Jval[r+1] <- Jcomputation$J
    exph.work <- Jcomputation$exp.signal
    
    ht.est.new <- h.eval(gamout,eval.pts,vshift)
    abs.diff <- abs(ht.est.new-ht.est)
    dh.abs.mean[r+1] <- mean(abs.diff)
    dh.abs.max[r+1] <- max(abs.diff)
    dh.abs.re[r+1] <- sum(abs.diff)/sum(abs(ht.est))
    ht.est <- ht.est.new
    
    #check: break
    if(identical(Jval[r+1],NaN)){
      stop(paste("NaN observed likelihood occurs at iteration: ",r))
    }
    Jascend[r+1] <- Jval[r+1]-Jval[r]
    if(Jascend[r+1]<thres && dh.abs.re[r+1]<thres){break}
  }
  
  ###calculate EL tuning
  ELtune <- data.frame(prop=lam.prop)
  exph <- exp(gamout$linear.predictors-vshift)
  eq.lam <- function(lambda){
    mean((exph-1)/(1+(exph-1)*lambda))
  }
  window <- (1/(n+m)-1)/(exph-1)
  window.l <- max(window[which(window<0)])
  window.u <- min(window[which(window>0)])
  ELtune$lam.root <- uniroot(Vectorize(eq.lam),c(window.l,window.u))$root
  
  ###compute F0 & F1:
  G0weight <- 1/(1+ELtune$prop*(exph-1))/(n+m)
  weight.dis <- data.frame(G0=G0weight,G1=exph*G0weight)#this one is approximately Ward's
  G0weight <- 1/(1+ELtune$lam.root*(exph-1))/(n+m)
  weight.dis.est <- data.frame(G0=G0weight,G1=exph*G0weight)#this one is more closed to CDF weights
  
  output <- list()
  output$fin.all <- gamout
  output$obj.df <- data.frame(J=Jval,J.ascend=Jascend,conv=conv,
                              dh.abs.re=dh.abs.re,dh.abs.mean=dh.abs.mean,dh.abs.max=dh.abs.max)
  output$conv <- ifelse(r<maxit,T,F)
  output$knots <- unique(knots.his)
  output$vshift <- vshift
  output$ELtune <- ELtune
  output$weight.Dis <- weight.dis
  output$weight.Dis.est <- weight.dis.est
  output$obsLL <- Jcomputation$obsLL
  output$nu <- nu
  output$exph <- exph#could be deleted
  output$G.l <- G.l#could be deleted
  return(output)
}
###calculate Q value (penalized log-likelihood multiplied by n+m)
obs.LL.pen <- function(id0,id1,gamout,nu,
                       pi0,pi1,lam.prop,vshift){
  ##################
  ###INPUTS:
  #id0,id1: subject indices for R=0 and R=1 
  #gamout: output of function gam
  #nu: smoothness penalty tuning parameter
  #pi0,pi1,lam.prop,vshift: some quantities required for observed loglik
  ###OUTPUTS: list of
  #J: value of observed loglik - penalty
  #exp.signal: exp(h(T_i)) for each i
  #obsLL: observed log-likelihood
  ##################
  exph <- exp(gamout$linear.predictors-vshift)
  obs.ll <- sum(log(pi0+(1-pi0)*exph[id0]))+sum(log(1-pi1+pi1*exph[id1]))-
    sum(log(1-lam.prop+lam.prop*exph))#compute l_*(h)
  b <- gamout$coefficients[-1]
  quadPen <- b%*%(gamout$smooth[[1]]$S[[1]])%*%b
  return(list(J=obs.ll-nu*quadPen,exp.signal=exph,obsLL=obs.ll))
}
###calculate h evaluated at dense grid
h.eval <- function(gamout,points,vshift){
  as.numeric(predict(gamout,newdata = data.frame(cov=points))-vshift)
}
###summarize some statistics from DRest0.EM
sum.DRest0.EM <- function(dat,EMout,s=0.2,ELtuneType=0,
                          eval.pts=NULL,s.range=c(0.1,0.3),
                          alt.auc=F){
  ##################
  ###INPUTS:
  #dat: dataframe consists of R and biomarker  
  #EMout: output of 'DRest0.EM'
  #s: where ROC is evaluated
  #ELtuneType: lambda in empirical likelihood: 0 for proportion; 1 for estimated root
  #s.range: where pAUC is computed
  #alt.auc: F: original method ; T: alternative method
  ###OUTPUTS:
  #stats.df: ROC,AUC,Youden; IC.df: aic bic; h.df: h(t) evaluations
  ##################
  ###ROC AUC:
  require(spatstat.univar)
  if(ELtuneType==0){
    Fg0 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis$G0)
    Fg1 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis$G1)
  }else{
    Fg0 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis.est$G0)
    Fg1 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis.est$G1)
  }
  auc <- (stieltjes(function(c){1-Fg1(c)},Fg0)[[1]]+stieltjes(function(c){Fg0(c)},Fg1)[[1]])/2
  stats.df <- data.frame(ROC=1-Fg1(quantile(Fg0,1-s)),AUC=auc)
  ###YI,cutoff
  require(rootSolve)
  rts <- uniroot.all(h.eval,c(min(dat$biomarker),max(dat$biomarker)),
                     gamout=EMout$fin.all,vshift=EMout$vshift)
  if(length(rts)==0){
    stats.df$cutoff <- NA;stats.df$Youden <- NA
  }else{
    YI.tem <- Fg0(rts)-Fg1(rts)
    idd <- which.max(YI.tem)
    stats.df$cutoff <- rts[idd]
    stats.df$Youden <- YI.tem[idd]
  }
  ###IC
  edf <- sum(EMout$fin.all$edf)
  edf1 <- sum(EMout$fin.all$edf1)
  IC.df <- data.frame(aic=2*(edf-EMout$obsLL),
                      bic=edf*log(nrow(dat))-2*(EMout$obsLL),
                      aic1=2*(edf1-EMout$obsLL),
                      bic1=edf1*log(nrow(dat))-2*(EMout$obsLL))
  ###h evaluation
  if(is.null(eval.pts)){
    eval.interval <- quantile(dat$biomarker,probs=c(0.05,0.95))
    eval.pts <- seq(from=eval.interval[1],to=eval.interval[2],length.out=30)
  }
  h.df <- data.frame(t=eval.pts,ht=h.eval(EMout$fin.all,eval.pts,EMout$vshift))
  
  ###pAUC:
  ROC.fun <- function(s){1-Fg1(quantile(Fg0,1-s))}
  res1 <- try((integrate(ROC.fun,s.range[1],s.range[2],
                         subdivisions = 500)$value)/(s.range[2]-s.range[1]),
              silent = T)
  if(class(res1)=="try-error"){
    x.eval <- seq(from=s.range[1],to=s.range[2],length.out=1000)
    res1 <- trapz(x.eval,ROC.fun(x.eval))/(s.range[2]-s.range[1])
  }
  stats.df$pAUC <- res1
  # b <- quantile(Fg0,1-s.range[1])
  # a <- quantile(Fg0,1-s.range[2])
  # res2 <- stieltjes(function(c){(1-Fg1(c))*(c>=a)*(c<=b)},Fg0)[[1]]/(s.range[2]-s.range[1])
  # res3 <- ROC.fun(s.range[1])*(1-s.range[1])-ROC.fun(s.range[2])*(1-s.range[2])+
  #   stieltjes(function(c){Fg0(c)*(c>=a)*(c<=b)},Fg1)[[1]]
  # res3 <- res3/(s.range[2]-s.range[1])
  # res4 <- (res2+res3)/2
  
  ###alternative auc:
  if(isTRUE(alt.auc)){
    res1 <- try((integrate(ROC.fun,0,1,subdivisions = 500)$value),silent = T)
    if(class(res1)=="try-error"){
      x.eval <- seq(from=0,to=1,length.out=1000)
      res1 <- trapz(x.eval,ROC.fun(x.eval))
    }
    stats.df$AUC <- res1
  }
  
 return(list(stats.df=stats.df,IC.df=IC.df,h.df=h.df#,pAUC=data.frame(d=res1,id1=res2,id2=res3,mid=res4)
             ))
}
est.ROC <- function(s.seq,dat,EMout,ELtuneType=0){
  require(spatstat.univar)
  if(ELtuneType==0){
    Fg0 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis$G0)
    Fg1 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis$G1)
  }else{
    Fg0 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis.est$G0)
    Fg1 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis.est$G1)
  }
  return(1-Fg1(quantile(Fg0,1-s.seq)))
}
###generate knots according to data quantiles.
gen.q.knots <- function(x,k,ord){
  #k: number of B-spline basis
  #ord: order of B-spline
  #generate inner knots:
  kts <- quantile(x,probs = seq(0,1,length.out=k-ord+1),names=F)
  #generate outer knots:
  dif <- diff(range(x))/(k-ord)
  for(jj in 1:ord){
    kts <- c(kts[1]-dif,kts,kts[length(kts)]+dif)
  }
  return(kts)
}
trapz <- function(x, y){
  idx <- 2:length(x)
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}






###for nu (tuning parameter) selection

##information criteria
IC.DRest.EM <- function(dat,pi0,pi1,
                        nu=1,k=30,ord=4,ord.pen=2,
                        maxit=500,thres=1e-4,
                        knots=NULL,ini.G=NULL){
  ##################
  ###INPUTS:
  #dat: dataframe consists of R and biomarker 
  #pi0=P(G=0|R=0); pi1=P(G=1|R=1)
  #nu: a seires of smoothness penalty tuning parameter
  #k: number of B-spline basis
  #ord, ord.pen: the orders of B-splines and the penalty
  #maxit: maximum number of iterations
  #thres: threshold value for stopping
  #knots: user specified knots.
  #ini.G: two components specifying initial pseudo-responses for Groups R=0,1
  
  ###OUTPUTS: a list of:
  #all.list: outputs of "DRest0.EM" for each tuning parameter
  #nu.selection: a dataframe recording tuning selection for different criteria
  #nu.res: the candidate tuning parameters and some outputs (IC values and convergence) 
  #################
  Nonu <- length(nu)
  all.list <- list()
  nu.res <- data.frame(tune=nu,conv=vector(length = Nonu),
                       aic=rep(NA,Nonu), bic=rep(NA,Nonu),
                       aic1=rep(NA,Nonu), bic1=rep(NA,Nonu),
                       edf=rep(NA,Nonu),edf1=rep(NA,Nonu))#,
                       #gcv=rep(NA,Nonu),gcv1=rep(NA,Nonu))
  
  # test0 <- try(DRest0.EM(dat,pi0,pi1,
  #                        nu=0.01,k=k,ord=ord,ord.pen=ord.pen,
  #                        maxit=maxit,thres=thres,knots=knots,ini.G=ini.G),
  #              silent = T)
  # if(!class(test0)=="try-error"){
  #   D0 <- test0$obsLL
  # }else{
  #   D0 <- NA
  # }
  
  
  for(kk in 1:Nonu){
    test <- try(DRest0.EM(dat,pi0,pi1,
                          nu=nu[kk],k=k,ord=ord,ord.pen=ord.pen,
                          maxit=maxit,thres=thres,knots=knots,ini.G=ini.G),
                silent = T)
    if(!class(test)=="try-error"){
      nu.res$conv[kk] <- test$conv
      ###IC
      edf <- sum(test$fin.all$edf)
      edf1 <- sum(test$fin.all$edf1)
      nu.res$edf[kk] <- edf;nu.res$edf1[kk] <- edf1
      nu.res$aic[kk] <- 2*(edf-test$obsLL)
      nu.res$aic1[kk] <- 2*(edf1-test$obsLL)
      nu.res$bic[kk] <- edf*log(nrow(dat))-2*(test$obsLL)
      nu.res$bic1[kk] <- edf1*log(nrow(dat))-2*(test$obsLL)
      
      all.list[[kk]] <- test
      
      #GCV
      # todiv <- (1-edf/(nrow(dat)))^2
      # nu.res$gcv[kk] <- (D0-test$obsLL)/todiv
      # todiv1 <- (1-edf1/(nrow(dat)))^2
      # nu.res$gcv1[kk] <- (D0-test$obsLL)/todiv1
    }
  }
  ###selection:
  nu.selection <- matrix(nrow=2,ncol=4)
  colnames(nu.selection) <- c("aic","bic","aic1","bic1")#,"gcv","gcv1")
  rownames(nu.selection) <- c("id","tune")
  for(jj in 1:4){
    id <- which.min(nu.res[,jj+2])
    if(length(id)>0){nu.selection[,jj] <- c(id,nu[id])}
  }
  return(list(all.list=all.list,tune.res=nu.res,
              tune.selection=as.data.frame(nu.selection)))
}


cv.DRest.EM <- function(dat,pi0,pi1,nfolds=5,seed=123,
                        nu=1,k=30,ord=4,ord.pen=2,
                        maxit=500,thres=1e-4,
                        knots=NULL,ini.G=NULL){
  ##################
  ###INPUTS:
  #dat: dataframe consists of R and biomarker 
  #pi0=P(G=0|R=0); pi1=P(G=1|R=1)
  #nfolds: number of folds for CV
  #seed: seed for random 
  #nu: a seires of smoothness penalty tuning parameter
  #k: number of B-spline basis
  #ord, ord.pen: the orders of B-splines and the penalty
  #maxit: maximum number of iterations
  #thres: threshold value for stopping
  #ini.G: two components specifying initial pseudo-responses for Groups R=0,1
  
  ###OUTPUTS: a list of:
  #nu.selection: a dataframe recording tuning selection for different criteria
  #nu.res: the candidate tuning parameters and some outputs (IC values and convergence) 
  #################
  id0 <- which(dat$R==0);id1 <- which(dat$R==1)
  n <- length(id0);m <- length(id1)
  lam.prop <- (pi1*m+(1-pi0)*n)/(n+m)#this quantity should be almost the same for the whole and train/test sample
  Nonu <- length(nu)
  grp0 <- rep_len(c(1:nfolds), length.out=n);grp1 <- rep_len(c(1:nfolds), length.out=m)
  
  set.seed(seed)
  cv.df0 <- data.frame(which=grp0,id=sample(c(id0),n))
  cv.df1 <- data.frame(which=grp1,id=sample(c(id1),m))
  dat$grp <- rbind(cv.df0[order(cv.df0$id),],cv.df1[order(cv.df1$id),])$which
  
  cv.oll <- cv.conv <- matrix(nrow=Nonu,ncol=nfolds)
  for(j in 1:nfolds){
    dat.test <- dat[dat$grp==j,]
    dat.train <- dat[dat$grp!=j,]
    if(length(ini.G)>2){
      ini.G.work <- ini.G[which(dat$grp!=j)]
    }else{ini.G.work <- ini.G}
    
    for(kk in 1:Nonu){
      test <- try(DRest0.EM(dat.train,pi0,pi1,
                            nu=nu[kk],k=k,ord=ord,ord.pen=ord.pen,
                            maxit=maxit,thres=thres,knots=knots,ini.G=ini.G.work),
                  silent = T)
      if(class(test)!="try-error"){
        #predict observed likelihood for test data!
        h.pred <- predict(test$fin.all,newdata = data.frame(cov=dat.test$biomarker))-test$vshift
        exph <- exp(as.numeric(h.pred))
        id0.test <- which(dat.test$R==0);id1.test <- which(dat.test$R==1)
        cv.oll[kk,j] <- sum(log(pi0+(1-pi0)*exph[id0.test]))+sum(log(1-pi1+pi1*exph[id1.test]))-
          sum(log(1-lam.prop+lam.prop*exph))#l_*(h) for test sample
        cv.conv[kk,j] <- test$conv
      }else{
        cv.conv[kk,j] <- F
      }
    }
  }
  conv.status <- apply(cv.conv, 1, mean)
  id.na <- which(conv.status!=1)
  oll <- apply(cv.oll,1, sum,na.rm=T)
  if(length(id.na)>0){oll[id.na] <- NA}
  df <- data.frame(tune=nu,oll.est=oll,conv=ifelse(conv.status==1,T,F))
  ###tuning selection:
  tune.id <- which.max(oll)
  return(list(tune=nu[tune.id],tune.id=tune.id,
              res=df,data.split=dat))
}























###EM Algorithm (parametric h)
DRest0.EM.para <- function(dat,pi0,pi1,degree=2,log1=F,log2=F,
                           maxit=500,thres=1e-4,ini.G=NULL){
  ##################
  ###INPUTS:
  #dat: dataframe consists of R and biomarker 
  #pi0=P(G=0|R=0); pi1=P(G=1|R=1)
  #degree: degree of the polynomial
  #log1,log2: if log(t), log(1-t) should be added as basis
  #maxit: maximum number of iterations
  #thres: threshold value for stopping
  #ini.G: two components specifying initial pseudo-responses for Groups R=0,1
  
  ###OUTPUTS:
  ##################
  cov <- dat$biomarker
  id0 <- which(dat$R==0);id1 <- which(dat$R==1)
  n <- length(id0);m <- length(id1)
  vshift <- log(pi1*m+(1-pi0)*n)-log((1-pi1)*m+pi0*n)
  lam.prop <- (pi1*m+(1-pi0)*n)/(n+m)
  conv <- vector(mode = "logical", length =0)
  Jval <-  Jascend <- numeric(length = 0)
  dh.abs.mean <- dh.abs.max <- dh.abs.re <- numeric(length = 0)
  ##design matrix:
  cov.mat <- matrix(1,nrow=n+m)
  if(degree>0){
    for(j in 1:degree){
      cov.mat <- cbind(cov.mat,cov^j)
    }
  }
  if(isTRUE(log1)){cov.mat <- cbind(cov.mat,log(cov))}
  if(isTRUE(log2)){cov.mat <- cbind(cov.mat,log(1-cov))}
  
  ###initialization:
  if(is.null(ini.G)){
    G.work <- ifelse(dat$R==1,pi1,1-pi0) #initial value for unknown G_i
  }else if(length(ini.G)==2){
    G.work <- ifelse(dat$R==1,ini.G[2],ini.G[1])
  }else{
    G.work <- ini.G
  }
  glmout0 <- glm(G.work~cov.mat-1,
                 family = binomial)#initial h coef estimation
  conv[1] <- glmout0$converged
  #CAUSION: I use Ghat as response directly; the weights are constant
  #compute initial J value (multiplied by n+m)
  Jcomputation <- obs.LL.pen.para(id0,id1,glmout0,pi0,pi1,lam.prop,vshift)
  Jval[1] <- Jcomputation$J
  exph.work <- Jcomputation$exp.signal
  #h evaluations:
  eval.interval <- quantile(dat$biomarker,probs=c(0.05,0.95))
  eval.pts <- seq(from=eval.interval[1],to=eval.interval[2],length.out=30)
  ht.est <- h.eval.para(glmout0,eval.pts,
                        degree=degree,log1=log1,log2=log2,vshift)
  
  ###EM algorithm
  for(r in 1:maxit){
    #E-step
    G.work[id0] <- (1-pi0)*exph.work[id0]/((1-pi0)*exph.work[id0]+pi0)
    G.work[id1] <- pi1*exph.work[id1]/(pi1*exph.work[id1]+1-pi1)
    G.work <- pmax(G.work,0)
    G.work <- pmin(G.work,1)
    glmout <- glm(G.work~cov.mat-1,family = binomial)
    conv[r+1] <- glmout$converged
    
    Jcomputation <- obs.LL.pen.para(id0,id1,glmout,pi0,pi1,lam.prop,vshift)
    Jval[r+1] <- Jcomputation$J
    exph.work <- Jcomputation$exp.signal
    
    ht.est.new <- h.eval.para(glmout,eval.pts,
                              degree=degree,log1=log1,log2=log2,vshift)
    abs.diff <- abs(ht.est.new-ht.est)
    dh.abs.mean[r+1] <- mean(abs.diff)
    dh.abs.max[r+1] <- max(abs.diff)
    dh.abs.re[r+1] <- sum(abs.diff)/sum(abs(ht.est))
    ht.est <- ht.est.new
    
    #check: break
    Jascend[r+1] <- Jval[r+1]-Jval[r]
    if(Jascend[r+1]<thres && dh.abs.re[r+1]<thres){break}
  }
  
  ###calculate EL tuning
  ELtune <- data.frame(prop=lam.prop)
  exph <- exp(glmout$linear.predictors-vshift)
  eq.lam <- function(lambda){
    mean((exph-1)/(1+(exph-1)*lambda))
  }
  window <- (1/(n+m)-1)/(exph-1)
  window.l <- max(window[which(window<0)])
  window.u <- min(window[which(window>0)])
  ELtune$lam.root <- uniroot(Vectorize(eq.lam),c(window.l,window.u))$root
  ###compute F0 & F1:
  G0weight <- 1/(1+ELtune$prop*(exph-1))/(n+m)
  weight.dis <- data.frame(G0=G0weight,G1=exph*G0weight)#this one is approximately Ward's
  G0weight <- 1/(1+ELtune$lam.root*(exph-1))/(n+m)
  weight.dis.est <- data.frame(G0=G0weight,G1=exph*G0weight)#this one is more closed to CDF weights
  
  output <- list()
  output$fin.all <- glmout
  output$obj.df <- data.frame(J=Jval,J.ascend=Jascend,conv=conv,
                              dh.abs.re=dh.abs.re,dh.abs.mean=dh.abs.mean,dh.abs.max=dh.abs.max)
  output$conv <- ifelse(r<maxit,T,F)
  output$vshift <- vshift
  output$ELtune <- ELtune
  output$weight.Dis <- weight.dis
  output$weight.Dis.est <- weight.dis.est
  output$obsLL <- Jcomputation$obsLL
  output$exph <- exph
  
  #final G update (serve as initials for semi-para method)
  G.work[id0] <- (1-pi0)*exph[id0]/((1-pi0)*exph[id0]+pi0)
  G.work[id1] <- pi1*exph[id1]/(pi1*exph[id1]+1-pi1)
  G.work <- pmax(G.work,0)
  G.work <- pmin(G.work,1)
  output$G.fin <- G.work
  return(output)
}
obs.LL.pen.para <- function(id0,id1,glmout,
                            pi0,pi1,lam.prop,vshift){
  ##################
  ###INPUTS:
  #id0,id1: subject indices for R=0 and R=1 
  #glmout: output of function glm
  #pi0,pi1,lam.prop,vshift: some quantities required for observed loglik
  ###OUTPUTS: list of
  #J: value of observed loglik 
  #exp.signal: exp(h(T_i)) for each i
  ##################
  exph <- exp(glmout$linear.predictors-vshift)
  obs.ll <- sum(log(pi0+(1-pi0)*exph[id0]))+sum(log(1-pi1+pi1*exph[id1]))-
    sum(log(1-lam.prop+lam.prop*exph))#compute l_*(h)
  return(list(J=obs.ll,exp.signal=exph))
}
h.eval.para <- function(glmout,eval.pts,degree,log1,log2,vshift){
  ##design matrix:
  cov.mat <- matrix(1,nrow=length(eval.pts))
  if(degree>0){
    for(j in 1:degree){
      cov.mat <- cbind(cov.mat,eval.pts^j)
    }
  }
  if(isTRUE(log1)){cov.mat <- cbind(cov.mat,log(eval.pts))}
  if(isTRUE(log2)){cov.mat <- cbind(cov.mat,log(1-eval.pts))}
  return(as.numeric(cov.mat%*%(glmout$coefficients))-vshift)
}
sum.DRest0para.EM <- function(dat,EMout,degree=2,log1=F,log2=F,
                              s=0.2,ELtuneType=0,eval.pts=NULL){
  ##################
  ###INPUTS:
  #dat: dataframe consists of R and biomarker  
  #EMout: output of 'DRest0.EM.para'
  #s: where ROC is evaluated
  #ELtuneType: lambda in empirical likelihood: 0 for proportion; 1 for estimated root
  ###OUTPUTS:
  #stats.df: ROC,AUC,Youden;  h.df: h(t) evaluations
  ##################
  ###ROC AUC:
  require(spatstat.univar)
  if(ELtuneType==0){
    Fg0 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis$G0)
    Fg1 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis$G1)
  }else{
    Fg0 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis.est$G0)
    Fg1 <- ewcdf(dat$biomarker,weights = EMout$weight.Dis.est$G1)
  }
  auc <- (stieltjes(function(c){1-Fg1(c)},Fg0)[[1]]+stieltjes(function(c){Fg0(c)},Fg1)[[1]])/2
  stats.df <- data.frame(ROC=1-Fg1(quantile(Fg0,1-s)),
                         AUC=auc)
  
  ###YI,cutoff
  require(rootSolve)
  rts <- uniroot.all(h.eval.para,c(min(dat$biomarker),max(dat$biomarker)),
                     glmout=EMout$fin.all,
                     degree=degree,log1=log1,log2=log2,vshift=EMout$vshift)
  if(length(rts)==0){
    stats.df$cutoff <- NA;stats.df$Youden <- NA
  }else{
    YI.tem <- Fg0(rts)-Fg1(rts)
    idd <- which.max(YI.tem)
    stats.df$cutoff <- rts[idd]
    stats.df$Youden <- YI.tem[idd]
  }
  
  ###h evaluation
  if(is.null(eval.pts)){
    eval.interval <- quantile(dat$biomarker,probs=c(0.05,0.95))
    eval.pts <- seq(from=eval.interval[1],to=eval.interval[2],length.out=30)
  }
  h.df <- data.frame(t=eval.pts,
                     ht=h.eval.para(EMout$fin.all,eval.pts,
                                    degree=degree,log1=log1,log2=log2,
                                    vshift=EMout$vshift))
  
  return(list(stats.df=stats.df,h.df=h.df))
}


get.init.G <- function(dat,pi0,pi1,degree=1,log1=F,log2=F,
                       maxit=500,thres=1e-5,ini.G=NULL){
  ##################
  ###INPUTS:
  #dat: dataframe consists of R and biomarker(single)
  #pi0=P(G=0|R=0); pi1=P(G=1|R=1)
  #degree: degree of the polynomial
  #log1,log2: if log(t), log(1-t) should be added as basis
  #maxit: maximum number of iterations
  #thres: threshold value for stopping
  #ini.G: two components specifying initial pseudo-responses for Groups R=0,1
  
  ###INPUTS: init.G: can be used in function "DRest0.EM"
  ##################
  tes <- try(DRest0.EM.para(dat,pi0=pi0,pi1=pi1,degree=degree,log1=log1,log2=log2,
                            maxit=maxit,thres=thres,ini.G=ini.G),silent = T)
  if(class(tes)!="try-error"){
    return(list(init.G=tes$G.fin,conv=tes$conv,glmout=tes$fin.all))
  }else{
    return(list(init.G=NULL,conv=F,glmout=NULL))
  }
}
