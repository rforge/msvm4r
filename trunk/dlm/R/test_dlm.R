setwd("~/Documents/Projects/StockMarket/R/study1/")
load("dmcpl.RData")

dat <- subset(dmcpl,id==1)
y <- cbind(dat$dy,dat$r)
x <- as.matrix(dat[,c("c1","c2")])
nx <- ncol(x)

  xn <- array(0.0,dim=c(ny,ny*nx,nt))
  for(i in 1:nt) {
    for(j in 1:ny) {
      jid <- (j-1)*nx+1
      xn[j,jid:(jid+nx-1),i] <- x[i,]
    }
  }
  x <- xn
  nx <- dim(x)[2]
  nt <- dim(x)[3]
  A <- diag(nx)
  Q <- .005*diag(nx)
  Sigma <- diag(nx)
  ws <- rep(0.0,nx)
  R <- 30*diag(ny)
  lt <- 1
  bt <- 1
  et <- nt
filt <- dlm.filter(y,x,A,Q,R,ws,Sigma,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt)
smth <- dlm.smoother(y,x,A,ws,Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt)
em <- dlm.em(smth,y,x,A,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt)
for(i in 1:100) {
  filt <- dlm.filter(y,x,A,em$Q,em$R,em$ws,em$Sigma,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt)
  cat(filt$like,"\n")
  smth <- dlm.smoother(y,x,A,em$ws,em$Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt)
  em <- dlm.em(smth,y,x,A,nt=dim(x)[3],nx=dim(x)[2],lt=1,bt=1,et=nt)
}

w <- filt$w
P <- filt$P
H <- filt$H
L <- filt$L


