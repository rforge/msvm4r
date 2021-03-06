
mvdlrm.filter <- function(y,x,A,Q,R,ws,Sigma,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt) {
  # y = t*m matrix
  # x = m*d*t array
  # A = d*d matrix
  # Q = d*d matrix
  # R = m*m matrix
  # ws = 
  # used by dlrm
  ny <- ncol(y)
  II <- diag(nx)
  w <- matrix(0.0,nrow=nt,ncol=nx)
  L <- P <- array(0.0,dim=c(nx,nx,nt))
  Var <- H <- array(0.0,dim=c(ny,ny,nt))
  onestep <- array(0.0,dim=c(nt,ny))
  K <- matrix(0.0,ncol=ny,nrow=nx)
  like <- 0.0
  for(case in 1:lt) {
    w[bt[case],] <- A%*%ws
    P[,,bt[case]] <- A%*%Sigma%*%t(A) + Q
    Var[,,bt[case]] <- x[,,bt[case]]%*%P[,,bt[case]]%*%t(x[,,bt[case]]) + R
    #if(!isSymmetric(Var[,,bt[case]])) Var[,,i] <- (t(Var[,,i])+Var[,,i])/2
    for(i in bt[case]:(et[case]-1)) {
      if(!isSymmetric(Var[,,i])) Var[,,i] <- (t(Var[,,i])+Var[,,i])/2
      H[,,i] <- solve(Var[,,i])
      K <- P[,,i]%*%t(x[,,i])%*%H[,,i]
      L[,,i] <- A%*%(II - K%*%x[,,i])
      onestep[i,] <- x[,,i]%*%as.matrix(w[i,])
      if(!any(is.na(y[i,]))) like <- like + dmvnorm(y[i,],onestep[i,],Var[,,i],log=TRUE)

      if(!any(is.na(y[i,]))) w[i+1,] <- A%*%K%*%y[i,] + L[,,i]%*%w[i,] else w[i+1,] <- A%*%w[i,]
      if(!any(is.na(y[i,]))) P[,,i+1] <- L[,,i]%*%P[,,i]%*%t(A) + Q else P[i+1,] <- A%*%P[,,i]%*%t(A) + Q
      Var[,,i+1] <- x[,,i+1]%*%P[,,i+1]%*%t(x[,,i+1]) + R
    }
    H[,,et[case]] <- solve(Var[,,et[case]])
    K <- P[,,et[case]]%*%t(x[,,et[case]])%*%H[,,et[case]]
    L[,,et[case]] <- A%*%(II - K%*%x[,,et[case]])
    
    onestep[et[case],] <- x[,,et[case]]%*%w[et[case],]
    if(!any(is.na(y[et[case],]))) like <- like + dmvnorm(y[et[case],],onestep[et[case],],Var[,,et[case]],log=TRUE)
  }
  return(list(w=w,P=P,H=H,L=L,onestep=onestep,like=like,var=Var))
}


mvdlrm.smoother <- function(y,x,A,ws,Sigma,w,P,H,L,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt) {
  # used by dlrm
  II <- diag(nx)
  wks <- matrix(0.0,nrow=nt,ncol=nx)
  Pks <- PCks <- array(0.0,dim=c(nx,nx,nt))
  P0s <- array(0.0,dim=c(nx,nx,lt))
  w0s <- array(0.0,dim=c(nx,lt))
  #t_P <- t_L <- matrix(0.0,nrow=nx,ncol=nx)  
  for(case in 1:lt) {
    u <- rep(0.0,nx)
    U <- matrix(0.0,ncol=nx,nrow=nx)
    for(i in et[case]:bt[case]) {
      #t_P <- matrix(P[i,],ncol=nx)
      #t_L <- matrix(L[i,],ncol=nx)
      if(!any(is.na(y[i,]))) u <- t(x[,,i])%*%H[,,i]%*%(y[i,] - x[,,i]%*%w[i,]) + t(L[,,i])%*%u else u <- t(A)%*%u
      if(!any(is.na(y[i,]))) U  <- t(x[,,i])%*%H[,,i]%*%x[,,i] + t(L[,,i])%*%U%*%L[,,i] else U <- t(A)%*%U%*%A
      wks[i,] <- w[i,] + P[,,i]%*%u
      Pks[,,i] <- P[,,i] - P[,,i]%*%U%*%P[,,i]
      ifelse(i>bt[case],PCks[,,i] <- (II - P[,,i]%*%U)%*%L[,,i-1]%*%P[,,i-1],PCks[,,i] <- (II - P[,,bt[case]]%*%U)%*%A%*%Sigma)
    }
    w0s[,case] <- ws + Sigma%*%t(A)%*%u
    P0s[,,case] <- Sigma - Sigma%*%t(A)%*%U%*%A%*%Sigma
  }
  return(list(wks=wks,Pks=Pks,PCks=PCks,w0s=w0s,P0s=P0s))
}

mvdlrm.em <- function(smth,y,x,A,nt=nrow(x),nx=ncol(x),lt=1,bt=1,et=nt) {
  ny <- ncol(y)
  a <- b <- c <- d <- matrix(0,ncol=nx,nrow=nx)
  for(case in 1:lt) {
    a <- a + smth$P0s[,,case] + smth$w0s[,case]%*%t(smth$w0s[,case])
    b <- b + smth$PCks[,,bt[case]] + smth$wks[bt[case],]%*%t(smth$w0s[,case])
    c <- c + smth$Pks[,,bt[case]] + smth$wks[bt[case],]%*%t(smth$wks[bt[case],])
    for(i in (bt[case]+1):et[case]) {
      a <- a + smth$Pks[,,i-1] + smth$wks[i-1,]%*%t(smth$wks[i-1,])
      b <- b + smth$PCks[,,i] + smth$wks[i,]%*%t(smth$wks[i-1,])
      c <- c + smth$Pks[,,i] + smth$wks[i,]%*%t(smth$wks[i,])
    }
  }
  
  p <- y - t(apply(x*array(rep(t(smth$wks),each=2),dim=c(ny,nx,nt)),c(1,3),sum))

  ws <- apply(smth$w0s,1,mean)

  for(case in 1:lt) {
    d <- d + (smth$w0s[,case] - ws)%*%t(smth$w0s[,case] - ws)
  }
  d <- d/lt
  Sigma <- apply(smth$P0s,c(1,2),mean) + d # see Max Welling, The Kalman Filter

  
  Q <- (1/nt)*(c - b%*%t(A) - A%*%t(b) + A%*%a%*%t(A))
  #Q <- (1/nt)*(c - tcrossprod(b,A) - tcrossprod(A,b) + tcrossprod(A%*%a),A))
  Q <- (t(Q)+Q)/2 # ensure symmetry
  
  A <- b%*%solve(a) # A = Phi
  if(sum(is.na(y)) == 0) {
    R <- 0
    #R <-(1/nt)*sum(p^2 + colSums(matrix((as.vector(apply(x,1,rep,times=nx))*as.vector(apply(x,1,rep,each=nx))),nrow=nx*nx)*matrix(smth$Pks,nrow=nx*nx)))
    #R <- (1/nt)* matrix(rowSums(apply(p,1,function(x) x%o%x)),ncol=2)
    R <- 0
    for(i in 1:nt) {
        R <- R + p[i,]%o%p[i,] + x[,,i]%*%smth$Pks[,,i]%*%t(x[,,i])
    }
    R <- R/nt
   

  } else {
    nna <- !is.na(y)
    nt <- sum(nna)
    p <- p[nna,]
    x <- x[nna,]
    smth$Pks <- smth$Pks[,,nna]
    #R <-(1/nt)*sum(p^2 + colSums(matrix((as.vector(apply(x,1,rep,times=nx))*as.vector(apply(x,1,rep,each=nx))),nrow=nx*nx)*matrix(smth$Pks,nrow=nx*nx)))
    R <- 0
    for(i in 1:nt) {
        R <- R + p[i,]%o%p[i,] + x[,,i]%*%smth$Pks[,,i]%*%t(x[,,1])
    }
    R <- R/nt
  }
  return(list(A=A,Q=Q,R=R,ws=ws,Sigma=Sigma))
}

mvdlrm.opt <- function(y,x,A,Q,R,ws,Sigma,Q.diag=FALSE,Sigma.diag=FALSE,R.diag=FALSE,est.ws=TRUE,est.Sigma=TRUE,est.A=TRUE,est.Q=TRUE,est.R=TRUE,method="BFGS",lt=1,bt=1,et=nt,hessian=FALSE,...) {
  # Dynamic Linear Regression Model
  #    using Kalman filter/smoother and EM

  # TODO: allow for general Q.c and Sigma.c

  dimR <- attr(R,"dim")
  dimA <- attr(A,"dim")
  dimws <- attr(ws,"dim")
  nx <- dim(x)[2]
  nt <- dim(x)[3]
  
  func <- function(par,lt,bt,et) {
    names <- names(par)
    if(length(tmp <- grep("A",names)) > 0) {
      A <- par[tmp]
      attr(A,"dim") <- dimA
    }
    if(length(tmp <- grep("Q",names)) > 0) {
      chk <- FALSE
      try({
        if(Q.diag) Q <- as.matrix(nlme::pdDiag(par[tmp])) else Q <- as.matrix(nlme::pdSymm(par[tmp]))
        chk <- TRUE
      })
      if(!chk) return(NA)
    }
    if(length(tmp <- grep("ws",names)) > 0) {
      ws <- par[tmp]
      attr(ws,"dim") <- dimws
    }
    if(length(tmp <- grep("R",names)) > 0) {
      chk <- FALSE
      try({
        if(R.diag) R <- as.matrix(nlme::pdDiag(par[tmp])) else R <- as.matrix(nlme::pdSymm(par[tmp]))
        chk <- TRUE
      })
      if(!chk) return(NA)
      #attr(R,"dim") <- dimR
      if(any(R==Inf)) return(NA)
    }
    if(length(tmp <- grep("Sigma",names)) > 0) {
      chk <- FALSE
      try({
        if(Sigma.diag) Sigma <- as.matrix(nlme::pdDiag(par[tmp])) else Sigma <- as.matrix(nlme::pdSymm(par[tmp]))
        chk <- TRUE
      })
      if(!chk) return(NA)
    }
    -mvdlrm.filter(y,x,A,Q,R,ws,Sigma,nt,nx,lt,bt,et)$like
  }

  start <- list()
  if(est.A) start$A <- as.numeric(A)
  if(est.Q) if(Q.diag) start$Q <- coef(nlme::pdDiag(Q)) else start$Q <- coef(nlme::pdSymm(Q))
  if(est.R) if(R.diag) start$R <- coef(nlme::pdDiag(R)) else start$R <- coef(nlme::pdSymm(R))
  if(est.ws) start$ws <- as.numeric(ws)
  if(est.Sigma) if(Sigma.diag) start$Sigma <- coef(nlme::pdDiag(Sigma)) else start$Sigma <- coef(nlme::pdSymm(Sigma))

  fit <- optim(unlist(start),func,method=method,hessian=hessian,lt=lt,bt=bt,et=et,...)

  names <- names(fit$par)
  if(length(tmp <- grep("A",names)) > 0) {
      A <- fit$par[tmp]
      attr(A,"dim") <- dimA
  }
  if(length(tmp <- grep("Q",names)) > 0) {
      Q <- fit$par[tmp]
      if(Q.diag) Q <- as.matrix(nlme::pdDiag(Q)) else Q <- as.matrix(nlme::pdSymm(Q))
  }
  if(length(tmp <- grep("ws",names)) > 0) {
      ws <- fit$par[tmp]
      attr(ws,"dim") <- dimws
  }
  if(length(tmp <- grep("R",names)) > 0) {
      R <- fit$par[tmp]
      if(R.diag) R <- as.matrix(nlme::pdDiag(R)) else R <- as.matrix(nlme::pdSymm(R))
  }
  if(length(tmp <- grep("Sigma",names)) > 0) {
      Sigma <- fit$par[tmp]
      if(Sigma.diag) Sigma <- as.matrix(nlme::pdDiag(Sigma)) else Sigma <- as.matrix(nlme::pdSymm(Sigma))
  }
  hessian <- fit$hessian
  return(list(A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,hessian=hessian))
}


mvdlrm <- function(y,x,maxit=100,ws,Sigma,A,Q,R,Q.c=NULL,Sigma.c=Q.c,R.c = NULL,ntimes=NULL,tol=1e-5,est=TRUE,est.ws=TRUE,est.Sigma=TRUE,est.A=TRUE,est.Q=TRUE,est.R=TRUE,filter.only=FALSE,verbose=FALSE,criterion=c("logLik","parameter"),method="BFGS",switch.LL=.5,switch.wait=5) {
  # Dynamic Linear Regression Model 
  #    using Kalman filter/smoother and EM/numerical optimization
  # author: M. Speekenbrink
  # version: 0.4
  # date: 26 Januari 2007
  # adapted from: Wu, L.S., Pai, J.S & Hosking, J.R.M. (1996). An algorithm for estimating parameters of state-space models. Statistics and Probability Letters, 28, 99-106.
  # (note different notation: x = M_t and A = Phi)
 
  criterion <- match.arg(criterion)
  if(length(criterion)!=1) stop("supplied value for `criterion' is unclear")
  
  if(is.vector(y)) y <- matrix(y,nrow=length(y))
  if(!is.matrix(y)) y <- as.matrix(y)
  ny <- ncol(y)
  y.names <- colnames(y)
  
  if(length(dim(x)) != 3) stop("x should be an array with dimension ny*nx*nt")
  x.names <- dimnames(x)[[2]]
  
  nx <- dim(x)[2]
  nt <- dim(x)[3]
  
  if(is.null(ntimes)) ntimes <- nt
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)

  if(missing(Sigma)) Sigma <- diag(nx)
  if(missing(Q)) Q <- diag(nx)
  if(missing(A)) A <- diag(nx)
  if(missing(R)) R <- diag(ny)
  if(missing(ws)) ws <- rep(0,nx)

  if(is.null(Q.c)) Q.c <- matrix(1,ncol=ncol(Q),nrow=ncol(Q))
  if(is.null(Sigma.c)) Sigma.c <- matrix(1,ncol=ncol(Sigma),nrow=nrow(Sigma))
  if(is.null(R.c)) R.c <- matrix(1,ncol=ncol(R),nrow=nrow(R))
  # ensure upper tri is identical to lower tri
  Q.c[upper.tri(Q.c)] <- t(Q.c)[upper.tri(Q.c)] 
  Sigma.c[upper.tri(Sigma.c)] <- t(Sigma.c)[upper.tri(Sigma.c)]
  R.c[upper.tri(R.c)] <- t(R.c)[upper.tri(R.c)]
  # check whether a diagonal is wanted
  if(sum(Q.c) == sum(diag(Q.c))) Q.diag <- TRUE else Q.diag <- FALSE
  if(sum(Sigma.c) == sum(diag(Sigma.c))) Sigma.diag <- TRUE else Sigma.diag <- FALSE
  if(sum(R.c) == sum(diag(R.c))) R.diag <- TRUE else R.diag <- FALSE

  filt <- mvdlrm.filter(y=y,x=x,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
  LL.old <- LL <- filt$like
  if(filter.only) {
    # skip smoother and break
    maxit <- -1
    smth <- list()
    smth$wks <- filt$w
    smth$Pks <- filt$P
  } else {
    smth <- mvdlrm.smoother(y=y,x=x,A=A,ws=ws,Sigma=Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
  }
  
  j <- 0
  k <- 0
  LL.dif <- tol+3
  
  if(est) {
    converge <- FALSE 
    if(verbose) cat("Kalman filter EM \n")
    opt.ok <- FALSE
    force.opt <- FALSE
  } else converge <- TRUE
  
  Hess <- NULL

  while(j <= maxit && !converge) {
    #abs(LL.dif) > tol
    # Expectation Maximisation
    
    j <- j+1

    if(((abs(LL.dif) < switch.LL) & opt.ok) || force.opt) {
      if(verbose) cat("starting optim \n")
      em <- mvdlrm.opt(y=y,x=x,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,Q.diag=Q.diag,Sigma.diag=Sigma.diag,R.diag=R.diag,est.ws=est.ws,est.Sigma=est.Sigma,est.A=est.A,est.Q=est.Q,est.R=est.R,method=method,hessian=hessian,lt=lt,bt=bt,et=et)
      opt.ok <- FALSE # avoid consecutive numerical optimization
      k <- 0
      converge <- TRUE # delete me
    } else {
      em <- mvdlrm.em(smth=smth,y=y,x=x,A=A,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
      #opt.ok <- TRUE #numerical optimization after em step ok
      k <- k+1
    }
    
    filt <- mvdlrm.filter(y=y,x=x,
      A = if(est.A) em$A else A,
      Q = if(est.Q) abs(1-Q.c)*Q + Q.c*em$Q else Q,
      R = if(est.R) abs(1-R.c)*R + R.c*em$R else R,
      ws = if(est.ws) em$ws else ws,
      Sigma = if(est.Sigma) abs(1-Sigma.c)*Sigma + Sigma.c*em$Sigma else Sigma,
      nt=nt,nx=nx,lt=lt,bt=bt,et=et)
  
    LL <- filt$like 
    LL.dif <- LL.old - LL
    
    if(criterion=="parameter") {
      if(all(abs(c(ws-em$ws,Sigma-(abs(1-Sigma.c)*Sigma + Sigma.c*em$Sigma),A-em$A,Q-(abs(1-Q.c)*Q + Q.c*em$Q),R-(abs(1-R.c)*R + R.c*em$R))) < tol )) converge <- TRUE
    } else {
      if(abs(LL.dif) < tol) converge <- TRUE
    }
    
    if(LL.dif < 0) {
      if(est.ws) ws <- em$ws
      if(est.Sigma) Sigma <- abs(1-Sigma.c)*Sigma + Sigma.c*em$Sigma
      if(est.A) A <- em$A # A = Phi
      if(est.Q) Q <- abs(1-Q.c)*Q + Q.c*em$Q
      if(est.R) R <- abs(1-R.c)*R + R.c*em$R
      Hess <- em$hessian
      LL.old <- LL
    } else {
      warning("Likelihood went down")
      if(j>switch.wait+1) {
        LL <- LL.old
        break
      } else {
        force.opt <- TRUE
#        opt.ok <- TRUE # avoid consecutive numerical optimization
#        LL.dif <- switch.LL-.1
      }
      #converge <- TRUE
    }
    
    smth <- mvdlrm.smoother(y=y,x=x,A=A,ws=ws,Sigma=Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=nt,nx=nx,lt=lt,bt=bt,et=et)

    if(k >= switch.wait) opt.ok <- TRUE
    if(verbose) cat("iteration",j,": LL =",LL,":LLdif =",LL.dif,"\n")
  }
  if(!filter.only) {
    filt <- mvdlrm.filter(y=y,x=x,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
    smth <- mvdlrm.smoother(y=y,x=x,A=A,ws=ws,Sigma=Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
  }
  # Number of free parameters (EM)
  npar <- 0
  
  if(est) {
      if(est.ws) npar <- npar + length(ws)
      if(est.Sigma) npar <- npar + sum(Sigma.c[lower.tri(Sigma.c,diag=TRUE)])
      if(est.A) npar <- npar + length(A)
      if(est.Q) npar <- npar + sum(Q.c[lower.tri(Q.c,diag=TRUE)])
      if(est.R) npar <- npar + sum(R.c[lower.tri(R.c,diag=TRUE)])
      
      if(maxit==0 | filter.only) npar <- 0 # nothing was actually estimated
  }
  colnames(filt$w) <- colnames(smth$wks) <- x.names
  return(list(call=call,response=y,predictor=x,weight=smth$wks,predicted=t(apply(x*array(rep(t(smth$wks),each=2),dim=c(ny,nx,nt)),c(1,3),sum)),weight.cov=smth$Pks,weight.filter=filt$w,predicted.onestep=filt$onestep,predicted.onestep.var=filt$var,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,LL=filt$like,npar=npar,niter=j,convergence=LL.dif,hessian=Hess))
}
