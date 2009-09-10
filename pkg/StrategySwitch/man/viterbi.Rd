\name{viterbi}
\alias{viterbi}
\title{Viterbi algorithm to determine Maximum A Posteriori strategy sequences}
\description{
The Viterbi algorithm is used to compute Maximum A Posteriori strategy sequences.
}
\usage{
viterbi(mod,A.id=NULL)
}
\arguments{
\item{mod}{An (estimated) strategy switch model}
\item{A.id}{a vector with integers indicating for each subject in mod which transition matrix (element in the list A) applies.}
}
\details{The Viterbi algorithm.
}
\value{A vector with integers indicating the MAP strategy for each trial.}
\author{Maarten Speekenbrink}
\references{
Speekenbrink, M., Lagnado, D. A, Wilkinson, L., Jahanshahi, M. & Shanks, D. R. (in press). Models of probabilistic category learning in Parkinson's disease: Strategy use and the effects of L-dopa. Journal of Mathematical Psychology. 
}
\examples{
## open weather prediction data
data(WPT)
## specify the Z matrix for the Constant Error version
Z <- rbind(rep(0,15),
            c(0,1,1,-1,-1,-2,0,0,0,0,2,1,1,-1,-1),
            c(-1,2,rep(0,13)),
            c(-1,0,2,rep(0,12)),
            c(1,0,0,-2,rep(0,11)),
            c(1,0,0,0,-2,rep(0,10)),
            c(0,1,1,-1,-1,-1,1,0,0,-1,1,0,2,-2,0))
## construct y and X lists
y <- X <- list()
for(i in 1:length(levels(WPT$id))) {
    dat <- subset(WPT,id==levels(WPT$id)[i])
    y[[i]] <- as.numeric(dat$r==0)
    X[[i]] <- model.matrix(r~x1*x2*x3*x4 - x1:x2:x3:x4,data=dat)
}
## set initial values for A and prior
A <- matrix(1/7,ncol=7,nrow=7)
prior <- c(1,rep(0,6))
## 
\dontrun{mod <- StrategySwitchModel(y=y,X=X,Z=Z,prior=prior,A=A,b=2.8,prior.est=rep(0,7),A.est=TRUE,b.est=TRUE,tol=1e-5,maxiter=2000,A.group=c(1,2),verbose=T,b.min=0.4054651)
## compute BIC
WPT$MAP <- viterbi(mod,A.id=c(1,2))}
}
