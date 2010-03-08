dlm.filter <- function(y,x,A,Q,R,ws,Sigma,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt) {
  # used by dlm
  # y: nt*ny matrix
  # x: ny*nx*nt array
  ny <- ncol(y)
  II <- diag(nx)
  w <- matrix(0.0,nrow=nt,ncol=nx)
  L <- P <- array(0.0,dim=c(ncol(Q),ncol(Q),nt))
  Cov <- H <- array(0.0,dim=c(ny,ny,nt))
  onestep <- matrix(0.0,nrow=nt,ncol=ny)
  K <- matrix(0.0,ncol=ny,nrow=nx)
  like <- 0.0
  for(case in 1:lt) {
    w[bt[case],] <- A%*%ws
    P[,,bt[case]] <- A%*%Sigma%*%t(A) + Q
    Cov[,,bt[case]] <- x[,,bt[case]]%*%P[,,bt[case]]%*%t(x[,,bt[case]]) + R
    for(i in bt[case]:(et[case]-1)) {
      H[,,i] <- solve(Cov[,,i])
      K <- P[,,i]%*%t(x[,,i])%*%H[,,i]
      L[,,i] <- A%*%(II - K%*%x[,,i])

      onestep[i,] <- x[,,i]%*%w[i,]
      if(!any(is.na(y[i,]))) like <- like + dmvnorm(y[i,],onestep[i,],Cov[,,i],log=TRUE)

      if(!any(is.na(y[i,]))) w[i+1,] <- A%*%K%*%y[i,] + L[,,i]%*%w[i,] else w[i+1,] <- A%*%w[i,]
      if(!any(is.na(y[i,]))) P[,,i+1] <- L[,,i]%*%P[,,i]%*%t(A) + Q else P[i+1,] <- A%*%P[,,i]%*%t(A) + Q
      Cov[,,i+1] <- x[,,i+1]%*%P[,,i+1]%*%t(x[,,i+1]) + R
    }
    H[,,et[case]] <- solve(Cov[et[case]])
    K <- P[,,et[case]]%*%t(x[,,et[case]])%*%H[,,et[case]]
    L[,,et[case]] <- A%*%(II - K%*%x[,,et[case]])
    
    onestep[et[case],] <- x[,,et[case]]%*%w[et[case],]
  }
  return(list(w=w,P=P,H=H,L=L,onestep=onestep,like=like,cov=Cov))
}

dlm.smoother <- function(y,x,A,ws,Sigma,w,P,H,L,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt) {
  # used by dlm
  II <- diag(nx)
  wks <- matrix(0.0,nrow=nt,ncol=nx)
  Pks <- PCks <- array(0.0,dim=c(nx,nx,nt))
  P0s <- array(0.0,dim=c(nx,nx,lt))
  w0s <- array(0.0,dim=c(nx,lt))
  for(case in 1:lt) {
    u <- rep(0.0,nx)
    U <- matrix(0.0,ncol=nx,nrow=nx)
    for(i in et[case]:bt[case]) {
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

dlm.em <- function(smth,y,x,A,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt) {
  a <- b <- c <- d <- matrix(0,ncol=nx,nrow=nx)
  for(case in 1:lt) {
    # used by dlrm
    a <- a + smth$P0s[,,case] + smth$w0s[,case]%*%t(smth$w0s[,case])
    b <- b + smth$PCks[,,bt[case]] + smth$wks[bt[case],]%*%t(smth$w0s[,case])
    c <- c + smth$Pks[,,bt[case]] + smth$wks[bt[case],]%*%t(smth$wks[bt[case],])
    for(i in (bt[case]+1):et[case]) {
      a <- a + smth$Pks[,,i-1] + smth$wks[i-1,]%*%t(smth$wks[i-1,])
      b <- b + smth$PCks[,,i] + smth$wks[i,]%*%t(smth$wks[i-1,])
      c <- c + smth$Pks[,,i] + smth$wks[i,]%*%t(smth$wks[i,])
    }
  }
  p <- y - apply(aperm(x,c(3,2,1))*as.vector(smth$wks),c(1,3),sum)
  ws <- apply(smth$w0s,1,mean)
  for(case in 1:lt) {
    d <- d + (smth$w0s[,case] - ws)%*%t(smth$w0s[,case] - ws)
  }
  d <- d/lt
  Sigma <- apply(smth$P0s,c(1,2),mean) + d # see Max Welling, The Kalman Filter

  Q <- (1/nt)*(c - b%*%t(A) - A%*%t(b) + A%*%a%*%t(A))
  Q <- (t(Q)+Q)/2 # ensure symmetry
  
  A <- b%*%solve(a) # A = Phi
  
  if(sum(is.na(y)) != 0) {
    nna <- !is.na(y)
    nt <- sum(nna)
    p <- p[nna,]
    x <- x[nna,]
    smth$Pks <- smth$Pks[,,nna]  
  }
  
  Rt <- matrix(0,0,ncol=ny,nrow=ny)
  for(i in 1:nt) {
    Rt <- Rt + p[i,]%*%t(p[i,]) + x[,,i]%*%smth$Pks[,,i]%*%t(x[,,i])
  }
  R <- (1/nt)*Rt
  return(list(A=A,Q=Q,R=R,ws=ws,Sigma=Sigma))
}

dlm <- function(formula,data,maxit=100,ws,Sigma,A,Q,R,Q.c=NULL,Sigma.c=Q.c,ntimes=NULL,tol=1e-5,est.ws=TRUE,est.Sigma=TRUE,est.A=TRUE,est.Q=TRUE,est.R=TRUE,filter.only=FALSE,verbose=FALSE,criterion=c("logLike","parameter"),method="BFGS",hessian=FALSE,switch.LL=.5,switch.wait=5) {
  # Dynamic Linear Regression Model 
  #    using Kalman filter/smoother and EM/numerical optimization
  # author: M. Speekenbrink
  # version: 0.4
  # date: 26 Januari 2007
  # adapted from: Wu, L.S., Pai, J.S & Hosking, J.R.M. (1996). An algorithm for estimating parameters of state-space models. Statistics and Probability Letters, 28, 99-106.
  # (note different notation: x = M_t and A = Phi)
 
  call <- match.call()
  criterion <- match.arg(criterion)
  if(length(criterion)!=1) stop("supplied value for `criterion' is unclear")
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf)

  if(is.vector(y)) y <- matrix(y,nrow=length(y))
  
  x <- if(!is.empty.model(mt)) {
    model.matrix(mt, mf, contrasts)
  } else matrix(, length(y), 0)
  x <-  matrix(as.numeric(x),nrow=dim(x)[1],ncol=dim(x)[2])
  
  ny <- ncol(y)
  nt <- nrow(x)
  nx <- ncol(x)
  
  # expand x
  xn <- array(0.0,dim=c(ny,ny*nx,nt))
  for(i in 1:nt) {
    for(j in 1:ny) {
      jid <- (j-1)*nx+1
      xn[j,jid:(jid+nx-1),i] <- x[i,]
    }
  }
  x <- xn
  
  #x.names <- dimnames(x)[[2]]
  #y.names <- levels(eval(attr(mt,"variables"),envir=mf)[[1]])
  
  nx <- dim(x)[2]
  
  if(is.null(ntimes)) ntimes <- nt
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)

  if(missing(Sigma)) Sigma <- diag(nx)
  if(missing(Q)) Q <- diag(nx)
  if(missing(A)) A <- diag(nx)
  if(missing(R)) R <- diag(ny)
  if(missing(ws)) ws <- rep(0,nx)

  if(is.null(Q.c)) Q.c <- matrix(1,ncol=ncol(Q),nrow=nrow(Q))
  if(is.null(Sigma.c)) Sigma.c <- matrix(1,ncol=ncol(Sigma),nrow=nrow(Sigma))
  # ensure upper tri is identical to lower tri
  Q.c[upper.tri(Q.c)] <- t(Q.c)[upper.tri(Q.c)] 
  Sigma.c[upper.tri(Sigma.c)] <- t(Sigma.c)[upper.tri(Sigma.c)]
  # check whether a diagonal is wanted
  if(sum(Q.c) == sum(diag(Q.c))) Q.diag <- TRUE else Q.diag <- FALSE
  if(sum(Sigma.c) == sum(diag(Sigma.c))) Sigma.diag <- TRUE else Sigma.diag <- FALSE

  filt <- dlm.filter(y=y,x=x,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
  LL.old <- LL <- filt$like
  if(filter.only) {
    # skip smoother and break
    maxit <- -1
    smth <- list()
    smth$wks <- filt$w
    smth$Pks <- filt$P
  } else {
    smth <- dlm.smoother(y=y,x=x,A=A,ws=ws,Sigma=Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
  }
  
  j <- 0
  k <- 0
  LL.dif <- tol+3
  
  if(verbose) cat("Kalman filter EM \n")
  converge <- FALSE
  opt.ok <- FALSE
  force.opt <- FALSE
  Hess <- NULL
  while(j <= maxit && !converge) {
    #abs(LL.dif) > tol
    # Expectation Maximisation
    
    j <- j+1

    if(((abs(LL.dif) < switch.LL) & opt.ok) || force.opt) {
      if(verbose) cat("starting optim \n")
      em <- dlm.opt(y=y,x=x,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,est.ws=est.ws,est.Sigma=est.Sigma,est.A=est.A,est.Q=est.Q,est.R=est.R,method=method,hessian=hessian,lt=lt,bt=bt,et=et)
      opt.ok <- FALSE # avoid consecutive numerical optimization
      k <- 0
      converge <- TRUE # delete me
    } else {
      em <- dlm.em(smth=smth,y=y,x=x,A=A,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
      #opt.ok <- TRUE #numerical optimization after em step ok
      k <- k+1
    }
    
    filt <- dlrm.filter(y=y,x=x,
      A = if(est.A) em$A else A,
      Q = if(est.Q) abs(1-Q.c)*Q + Q.c*em$Q else Q,
      R = if(est.R) em$R else R,
      ws = if(est.ws) em$ws else ws,
      Sigma = if(est.Sigma) abs(1-Sigma.c)*Sigma + Sigma.c*em$Sigma else Sigma,
      nt=nt,nx=nx,lt=lt,bt=bt,et=et)
  
    LL <- filt$like 
    LL.dif <- LL.old - LL
    
    if(criterion[1]==1) {
      if(all(abs(c(ws-em$ws,Sigma-(abs(1-Sigma.c)*Sigma + Sigma.c*em$Sigma),A-em$A,Q-(abs(1-Q.c)*Q + Q.c*em$Q),R-em$R)) < tol )) converge <- TRUE
    } else {
      if(abs(LL.dif) < tol) converge <- TRUE
    }
    
    if(LL.dif < 0) {
      if(est.ws) ws <- em$ws
      if(est.Sigma) Sigma <- abs(1-Sigma.c)*Sigma + Sigma.c*em$Sigma
      if(est.A) A <- em$A # A = Phi
      if(est.Q) Q <- abs(1-Q.c)*Q + Q.c*em$Q
      if(est.R) R <- em$R
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
    
    smth <- dlm.smoother(y=y,x=x,A=A,ws=ws,Sigma=Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=nt,nx=nx,lt=lt,bt=bt,et=et)

    if(k >= switch.wait) opt.ok <- TRUE
    if(verbose) cat("iteration",j,": LL =",LL,":LLdif =",LL.dif,"\n")
  }
  if(!filter.only) {
    filt <- dlm.filter(y=y,x=x,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
    smth <- dlm.smoother(y=y,x=x,A=A,ws=ws,Sigma=Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
  }
  # Number of free parameters (EM)
  npar <- 0
  if(est.ws) npar <- npar + length(ws)
  if(est.Sigma) npar <- npar + sum(Sigma.c[lower.tri(Sigma.c,diag=TRUE)])
  if(est.A) npar <- npar + length(A)
  if(est.Q) npar <- npar + sum(Q.c[lower.tri(Q.c,diag=TRUE)])
  if(est.R) npar <- npar + 1
  
  if(maxit==0 | filter.only) npar <- 0 # nothing was actually estimated
  
  colnames(filt$w) <- colnames(smth$wks) <- x.names
  return(list(call=call,response=y,predictor=x,weight=smth$wks,predicted=rowSums(smth$wks*x),weight.cov=smth$Pks,weight.filter=filt$w,predicted.onestep=filt$onestep,predicted.onestep.var=filt$var,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,LL=filt$like,npar=npar,niter=j,convergence=LL.dif,hessian=Hess))
}


dlm.opt <- function(y,x,A,Q,R,ws,Sigma,Q.diag=FALSE,Sigma.diag=FALSE,est.ws=TRUE,est.Sigma=TRUE,est.A=TRUE,est.Q=TRUE,est.R=TRUE,method="BFGS",lt=1,bt=1,et=nt,hessian=FALSE,...) {
  # Dynamic Linear Regression Model
  #    using Kalman filter/smoother and EM

  # TODO: allow for general Q.c and Sigma.c

  dimR <- attr(R,"dim")
  dimA <- attr(A,"dim")
  dimws <- attr(ws,"dim")
  nx <- ncol(x)
  nt <- nrow(x)
  
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
      R <- exp(par[tmp])
      attr(R,"dim") <- dimR
      if(R==Inf) return(NA)
    }
    if(length(tmp <- grep("Sigma",names)) > 0) {
      chk <- FALSE
      try({
        if(Sigma.diag) Sigma <- as.matrix(nlme::pdDiag(par[tmp])) else Sigma <- as.matrix(nlme::pdSymm(par[tmp]))
        chk <- TRUE
      })
      if(!chk) return(NA)
    }

    -dlrm.filter(y,x,A,Q,R,ws,Sigma,nt,nx,lt,bt,et)$like
  }

  start <- list()
  if(est.A) start$A <- as.numeric(A)
  if(est.Q) start$Q <- coef(nlme::pdSymm(Q))
  if(est.R) start$R <- coef(nlme::pdSymm(R))
  if(est.ws) start$ws <- as.numeric(ws)
  if(est.Sigma) start$Sigma <- coef(nlme::pdSymm(Sigma))

  fit <- optim(unlist(start),func,method=method,hessian=hessian,lt=lt,bt=bt,et=et,...)

  names <- names(fit$par)
  if(length(tmp <- grep("A",names)) > 0) {
      A <- fit$par[tmp]
      attr(A,"dim") <- dimA
  }
  if(length(tmp <- grep("Q",names)) > 0) {
      Q <- fit$par[tmp]
      Q <- as.matrix(nlme::pdSymm(Q))
  }
  if(length(tmp <- grep("ws",names)) > 0) {
      ws <- fit$par[tmp]
      attr(ws,"dim") <- dimws
  }
  if(length(tmp <- grep("R",names)) > 0) {
      R <- fit$par[tmp]
      R <- as.matrix(nlme::pdSymm(R))
  }
  if(length(tmp <- grep("Sigma",names)) > 0) {
      Sigma <- fit$par[tmp]
      Sigma <- as.matrix(nlme::pdSymm(Sigma))
  }
  hessian <- fit$hessian
  return(list(A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,hessian=hessian))
}

dlm.hess <- function(y,x,A,Q,R,ws,Sigma,Q.diag,Sigma.diag,est.ws,est.Sigma,est.A,est.Q,est.R,lt=1,bt=1,et=nt,method="numDeriv",...) {
  #require(numDeriv)
  if(method=="numDeriv") require(numDeriv)
  if(method=="fdHess") require(nlme)
  
  dimR <- attr(R,"dim")
  dimA <- attr(A,"dim")
  dimws <- attr(ws,"dim")
  nx <- ncol(x)
  nt <- nrow(x)

  func <- function(par,lt,bt,et) {
    names <- names(par)
    if(length(tmp <- grep("A",names)) > 0) {
      A <- par[tmp]
      attr(A,"dim") <- dimA
    }
    if(length(tmp <- grep("Q",names)) > 0) {
      chk <- FALSE
      try({
        if(Q.diag) Q <- diag(par[tmp]) else {
          Q <- matrix(0,ncol=nx,nrow=nx)
          Q[upper.tri(Q,diag=TRUE)] <- par[tmp]
          Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
        }
        chk <- TRUE
      })
      if(!chk) return(NA)
    }
    if(length(tmp <- grep("ws",names)) > 0) {
      ws <- par[tmp]
      attr(ws,"dim") <- dimws
    }
    if(length(tmp <- grep("R",names)) > 0) {
      R <- par[tmp]
      attr(R,"dim") <- dimR
    }
    if(length(tmp <- grep("Sigma",names)) > 0) {
      chk <- FALSE
      try({
        if(Sigma.diag) Sigma <- diag(par[tmp]) else {
          Sigma <- matrix(0,ncol=nx,nrow=nx)
          Sigma[upper.tri(Sigma,diag=TRUE)] <- par[tmp]
          Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
        }
        chk <- TRUE
      })
      if(!chk) return(NA)
    }

    -dlrm.filter(y,x,A,Q,R,ws,Sigma,nt,nx,lt,bt,et)$like
  }

  start <- list()
  if(est.A) start$A <- as.numeric(A)
  if(est.Q) {
    if(Q.diag) start$Q <- diag(Q) else start$Q <- Q[upper.tri(Q,diag=TRUE)]
  }
  if(est.R) start$R <- R
  if(est.ws) start$ws <- as.numeric(ws)
  if(est.Sigma) {
    if(Sigma.diag) start$Sigma <- diag(Sigma) else start$Sigma <- Sigma[upper.tri(Sigma,diag=TRUE)]
  }
  if(method=="fdHess") {
    hess <- fdHess(pars=unlist(start),fun=func,lt=lt,bt=bt,et=et,...)
  }
  if(method=="numDeriv") {
    hess <- hessian(func=func,x=unlist(start),method="Richardson", method.args=list(d=.01,r=6),lt=lt,bt=bt,et=et,...)
    colnames(hess) <- rownames(hess) <- names(unlist(start))
  }
  hess
}

dlm.Hessian <- function(mod,method="numDeriv",ntimes,...) {
  y <- mod$response
  x <- mod$predictor
  A <- mod$A
  Q <- mod$Q
  R <- mod$R
  Sigma <- mod$Sigma
  ws <- mod$ws

  est.A <- eval(mod$call$est.A)
  if(is.null(est.A)) est.A <- TRUE
  est.Q <- eval(mod$call$est.Q)
  if(is.null(est.Q)) est.Q <- TRUE
  est.R <- eval(mod$call$est.R)
  if(is.null(est.R)) est.R <- TRUE
  est.Sigma <- eval(mod$call$est.Sigma)
  if(is.null(est.Sigma)) est.Sigma <- TRUE
  est.ws <- eval(mod$call$est.ws)
  if(is.null(est.ws)) est.ws <- TRUE

  Q.c <- eval(mod$call$Q.c)
  if(is.null(Q.c)) Q.diag <- FALSE else {
    Q.c[upper.tri(Q.c)] <- t(Q.c)[upper.tri(Q.c)]
    if(sum(Q.c) == sum(diag(Q.c))) Q.diag <- TRUE else Q.diag <- FALSE
  }
  Sigma.c <- eval(mod$call$Sigma.c)
  if(is.null(Sigma.c)) Sigma.diag <- FALSE else {
    Sigma.c[upper.tri(Sigma.c)] <- t(Sigma.c)[upper.tri(Sigma.c)]
    if(sum(Sigma.c) == sum(diag(Sigma.c))) Sigma.diag <- TRUE else Sigma.diag <- FALSE
  }

  if(missing(ntimes)) ntimes <- eval(mod$call$ntimes)
  if(is.null(ntimes)) ntimes <- nrow(mod$response)
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)

  hess <- dlrm.hess(y=y,x=x,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,Q.diag=Q.diag,Sigma.diag=Sigma.diag,est.ws=est.ws,est.Sigma=est.Sigma,est.A=est.A,est.Q=est.Q,est.R=est.R,lt=lt,bt=bt,et=et,method=method,...)
  hess
}

dlrm.CI <- function(mod,p=.95) {
  conf.level <- p
  crit <- qnorm((1 + conf.level)/2)
  A <- mod$A
  Q <- mod$Q
  R <- mod$R
  Sigma <- mod$Sigma
  ws <- mod$ws

  est.A <- eval(mod$call$est.A)
  if(is.null(est.A)) est.A <- TRUE
  est.Q <- eval(mod$call$est.Q)
  if(is.null(est.Q)) est.Q <- TRUE
  est.R <- eval(mod$call$est.R)
  if(is.null(est.R)) est.R <- TRUE
  est.Sigma <- eval(mod$call$est.Sigma)
  if(is.null(est.Sigma)) est.Sigma <- TRUE
  est.ws <- eval(mod$call$est.ws)
  if(is.null(est.ws)) est.ws <- TRUE

  Q.c <- eval(mod$call$Q.c)
  if(is.null(Q.c)) Q.diag <- FALSE else {
    Q.c[upper.tri(Q.c)] <- t(Q.c)[upper.tri(Q.c)]
    if(sum(Q.c) == sum(diag(Q.c))) Q.diag <- TRUE else Q.diag <- FALSE
  }
  Sigma.c <- eval(mod$call$Sigma.c)
  if(is.null(Sigma.c)) Sigma.diag <- FALSE else {
    Sigma.c[upper.tri(Sigma.c)] <- t(Sigma.c)[upper.tri(Sigma.c)]
    if(sum(Sigma.c) == sum(diag(Sigma.c))) Sigma.diag <- TRUE else Sigma.diag <- FALSE
  }
  start <- list()
  if(est.A) start$A <- as.numeric(A)
  if(est.Q) {
    if(Q.diag) start$Q <- diag(Q) else start$Q <- Q[upper.tri(Q,diag=TRUE)]
  }
  if(est.R) start$R <- R
  if(est.ws) start$ws <- as.numeric(ws)
  if(est.Sigma) {
    if(Sigma.diag) start$Sigma <- diag(Sigma) else start$Sigma <- Sigma[upper.tri(Sigma,diag=TRUE)]
  }

  theta.hat <- unlist(start)
  inv.fish <- solve(mod$hessian)
  lower <- theta.hat - crit * sqrt(diag(inv.fish))
  upper <- theta.hat + crit * sqrt(diag(inv.fish))
  pm <- crit * sqrt(diag(inv.fish))
  cbind(theta.hat,lower,upper,pm)
}

