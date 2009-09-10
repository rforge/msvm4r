viterbi <- function(mod,A.id=NULL) {
    ntimes <- unlist(lapply(mod$B,nrow))
    nt <- sum(ntimes)
    lt <- length(ntimes)
		et <- cumsum(ntimes)
		bt <- c(1,et[-lt]+1)
    if(is.list(mod$A)) {
      if(length(mod$A) > 1) {
        if(is.null(A.id)) stop("must give A.id!")
        if(length(A.id)!=lt) stop("A.id must have length N")
      } else {
        A.id <- rep(1,lt)
      }
    } else {
      A.id <- rep(1,lt)
    }
    if(length(unique(A.id)) == 1) {
      map <- unlist(lapply(mod$B,rule.viterbi,A=mod$A,prior=mod$prior))
    } else {
      map <- vector()
      for(i in 1:length(mod$B)) {
        map <- c(map,rule.viterbi(A=mod$A[[A.id[i]]],B=mod$B[[i]],prior=mod$prior))
      }
    }
    map
}

rule.fb <- function(A,B,prior) {

    # Forward-Backward algorithm (used in Baum-Welch)
    # Returns alpha, beta, and full data likelihood
    # A = K*K matrix with transition probabilities, from row to column
    # B = T*K matrix with elements ab_{ij} = P(y_i|s_j)
    # pi = K vector with prior probabilities
    
    # NOTE: to prevent underflow, alpha and beta are scaled, using c
    
    nt <- nrow(B)
    ns <- ncol(A)
    alpha <- matrix(ncol=ns,nrow=nt)
    beta <- matrix(ncol=ns,nrow=nt)
    c <- vector(length=nt)
    
    alpha[1,] <- prior*B[1,] # initialize
    c[1] <- 1/sum(alpha[1,])
    alpha[1,] <- alpha[1,]*c[1]
    for(i in 1:(nt-1)) {
        alpha[i+1,] <- (t(A)%*%alpha[i,])*B[i+1,]
        c[i+1] <- 1/sum(alpha[i+1,])
        alpha[i+1,] <- c[i+1]*alpha[i+1,] 
    }
    
    beta[nt,] <- 1*c[nt] # initialize
    for(i in (nt-1):1) {
        beta[i,] <- (A%*%(B[i+1,]*beta[i+1,]))*c[i]
    }
    
    gamma <- alpha*beta/c
    
    xi <- array(dim=c(nt-1,ns,ns))
    for(i in 1:(nt-1)) {
        xi[i,,] <- (alpha[i,]%*%t(B[i+1,]*beta[i+1,]))*A
    }
    
    like <- -sum(log(c))
    return(list(alpha=alpha,beta=beta,gamma=gamma,xi=xi,c=c,logLike=like))
}

rule.viterbi <- function(A,B,prior) {
    # returns the most likely state sequence
    nt <- nrow(B)
    ns <- ncol(A)
    delta <- psi <- matrix(nrow=nt,ncol=ns)
    state <- vector(length=nt)
    # initialization
    delta[1,] <- - (log(prior) + log(B[1,]))
    psi[1,] <- 0
    # recursion
    for(i in 2:nt) {
        for(j in 1:ns) {
            delta[i,j] <- min(delta[i-1,] - log(A[,j])) - log(B[i,j])
            k <- which.min(delta[i-1,] - log(A[,j]))
            if(length(k) == 0) k <- 0
            psi[i,j] <- k
        }
    }
    #trace maximum likely state
    state[nt] <- which.min(delta[nt,])
    for(i in (nt-1):1) {
        state[i] <- psi[i+1,state[i+1]]
    }
    return(state)
}

freepar.A <- function(A.est) {
  A.id <- unique(as.numeric(A.est[A.est!=0]))
  ns <- ncol(A.est)
  mat <- matrix(0,nrow=ns,ncol=length(A.id))
  for(i in 1:ns) {
    for(j in 1:length(A.id)) mat[i,j] <- sum(A.est[i,]==A.id[j])
  }
  mat <- mat[rowSums(mat)!=0,]
  length(A.id) - qr(mat,rep(1,nrow(mat)))$rank
}

StrategySwitch <- function(y,X,Z,prior,A,b,tol=1e-4,maxiter=200,A.est=TRUE,prior.est,b.est,A.group=rep(1,length(y)),verbose=FALSE,b.min=-Inf) {
  #
  # Baum-Welch (EM) algorithm for ML estimates of the Strategy Switch model
  # M. Speekenbrink
  #
  # y = list of lenght N with response variable
  # X = list of length N, each element a T_i*k design matrix (predictors)
  # Z = S*k proportional weight matrix,
  # prior = S vector with prior probabilities
  #
  # A.est: matrix of dim(A.est)=dim(A) with integers 0,1,... where 0 indicates fixed and 1,2,... indicate unique free parameters
  # prior.est: similar indicator vector as for A.est
  # b.est: N vector with integers 0,1,... where 0 indicates fixed and 1,2,... indicate unique free parameters
  #        defaults to (1,2,...,N)
  ni <- length(y)
  ns <- nrow(Z)
  ng <- length(unique(A.group))
  
  A.group <- as.numeric(factor(A.group,labels=1:ng))
  if(length(A.group)!=ni) stop("A.group must have length",ni)
  
  tmp <- list()
  if(!is.list(A)) {
    for(i in 1:ng) { tmp[[i]] <- A }
    A <- tmp
  }
  if(length(A)!=ng) stop("A must have length",ng)
  
  if(length(b) < ni) b <- c(b,rep(b[length(b)],ni-length(b)))

    if(missing(A.est)) {
        estA <- FALSE
        A.id <- numeric(0)
    } else {
      if(is.logical(A.est)) {
        estA <- A.est
        if(estA) A.est <- t(matrix(1:(ns^2),ncol=ns))
      } else {
        estA <- TRUE
        if(length(A.est)!=ns^2) stop("A.est should be an ns*ns square matrix")
      }
    }
    if(estA) {
        A.id <- unique(as.numeric(A.est[A.est!=0]))
        A.fix <- list()
        for(i in 1:ng) {
          A.fix[[i]] <- A[[i]][A.est==0]
        }
    }
    if(missing(prior.est)) {
        estprior <- FALSE
        prior.id <- numeric(0)
    } else {
      if(is.logical(prior.est)) {
        estprior <- prior.est
        if(estprior) prior.est <- 1:ns
      } else {
        estprior <- TRUE
        if(length(prior.est)!=ns) stop("prior.est should have length ns")
      }
    }
    if(estprior) {
        prior.id <- unique(as.numeric(prior.est[prior.est!=0]))
        prior.fix <- prior[prior.est==0]
    }
    if(missing(b.est)) {
        estb <- FALSE
        b.id <- numeric(0)
    } else {
      if(is.logical(b.est)) {
        estb <- b.est
        if(estb) b.est <- 1:length(y)
      } else {
        estb <- TRUE
        if(length(b.est)!=length(b)) stop("b.est must have length",length(b))
      }
    }
    if(estb) {
        b.id <- order(unique(as.numeric(b.est[b.est!=0])))
        b.fix <- b[b.est==0]
    }

  dat <- fbo <- B <- list()
  LL <- nt <- vector(length=ni)

  A_num <- array(dim=c(ni,nrow(A[[1]]),ncol(A[[1]])))
  A_denom <- matrix(nrow=ni,ncol=ncol(A[[1]]))
  
  priors <- matrix(nrow=ni,ncol=length(prior))
    
  for(i in 1:ni) {
    nt[i] <- nrow(X[[i]])
    dat[[i]] <- data.frame(cbind(y=y[[i]],x=as.numeric(X[[i]]%*%t(Z))))
    py <- matrix(1/(1+exp(-b[i]*dat[[i]]$x)),ncol=ns)
    B[[i]] <- (py^y[[i]])*(1-py)^(1-y[[i]])
    fbo[[i]] <- rule.fb(A=A[[A.group[i]]],B=B[[i]],prior=prior)
    LL[i] <- fbo[[i]]$logLike
  }

  
  LL.dif <- tol+1
  j <- 0
  
  while(j <= maxiter && LL.dif > tol) {
    LL.old <- LL
    j <- j+1
    
    # updates
    
    
    if(estprior) {
      for(i in 1:ni) priors[i,] <- fbo[[i]]$gamma[1,]
      prior <- colMeans(priors)
      for(k in prior.id) prior <- replace(prior,prior.est == k,mean(prior[prior.est == k]))
      if(length(prior.fix > 0)) prior <- replace(prior,prior.est == 0,prior.fix)
      prior <- prior/sum(prior)
    }
    if(estA) {
      for(i in 1:ni) {
        A_num[i,,] <- apply(fbo[[i]]$xi,c(2,3),sum)
        A_denom[i,] <- colSums(fbo[[i]]$gamma[1:(nt[i]-1),])
      }
      for(i in 1:ng) {
	if(sum(A.group==i)>1) A[[i]] <- apply(A_num[A.group==i,,],c(2,3),sum)/colSums(A_denom[A.group==i,]) else A[[i]] <- A_num[A.group==i,,]/A_denom[A.group==i,]
      }
      for(i in 1:ng) {
        tps <- matrix(rowMeans(matrix(unlist(lapply(fbo[A.group==i],function(x) colMeans(x$gamma))),nrow=ns)),nrow=ns,ncol=ns)
        for(k in A.id) {
          tmp <- tps[A.est == k]
          tmp <- tmp/sum(tmp)
          A[[i]] <- replace(A[[i]],A.est == k,sum(tmp*A[[i]][A.est == k]))
        }
        if(length(A.fix[[i]]) > 0) A[[i]] <- replace(A[[i]],A.est == 0,A.fix[[i]])
        A[[i]] <- A[[i]]/rowSums(A[[i]])
      }
    }
    if(estb) {
      for(k in b.id) {
        tdat <- data.frame()
        gamma <- numeric(0)
        for(i in (1:ni)[b.est==k]) {
          tdat <- rbind(tdat,dat[[i]][-(1:nt[i]),]) # delete random strategy
          gamma <- c(gamma,as.numeric(fbo[[i]]$gamma[,-1])) # delete random strategy
        }
        tdat$w <- gamma
        b <- replace(b,b.est==k,glm(y~x-1,data=tdat,family=binomial(),weights=tdat$w)$coefficients)
      }
      if(length(b.fix)>0) b <- replace(b,b.est==0,b.fix)
      b[b<b.min] <- b.min
    }

    for(i in 1:ni) {
      py <- matrix(1/(1+exp(-b[i]*dat[[i]]$x)),ncol=ns)
      B[[i]] <- (py^y[[i]])*(1-py)^(1-y[[i]])
      fbo[[i]] <- rule.fb(A=A[[A.group[i]]],B=B[[i]],prior=prior)
      LL[i] <- fbo[[i]]$logLike
    }

    LL.dif <- sum(LL)-sum(LL.old)
    if(verbose) cat(paste("iteration",j,": difference =",LL.dif,"\n"))
    if(LL.dif < 0) LL <- LL.old
  }
  npar <- 0
  if(estprior) {
    if(length(prior.id) > 0) npar <- npar + length(prior.id)-1
  }
  if(estb) npar <- npar + length(b.id)
  if(estA) {
    # FIX THIS
    npar <- npar + ng*freepar.A(A.est)
  } else {
    if(estA) npar <- npar + ng*ns*(ns-1)
  }
  #npar <- sum(c(length(A.id),length(prior.id),length(b.id)))
  return(list(A=A,B=B,prior=prior,b=b,LL=LL,df=npar,iterations=j-1,convergence=LL.dif))
}
