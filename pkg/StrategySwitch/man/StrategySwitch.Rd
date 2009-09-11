\name{StrategySwitch}
\alias{StrategySwitch}
\title{Strategy Switch model}
\description{
Estimates the Strategy Switch model (Speekenbrink et al., 2009).
}
\usage{
StrategySwitch(y,X,Z,prior,A,b,tol=1e-4,maxiter=200,
A.est=TRUE,prior.est,b.est,A.group=rep(1,length(y)),
verbose=FALSE,b.min=-Inf)
}
\arguments{
\item{y}{A list with each element being a vector with binary responses of a subject.}
\item{X}{A list with each element a matrix of binary indicators for presence X[i,j] = 1 or absence X[i,j] = 0 of cue j at trial i of a subject.}
\item{Z}{A matrix with coefficients for }
\item{prior}{A vector with (starting values) for the prior probabilties of the strategies.}
\item{A}{A list with starting values for the transition matrices.}
\item{b}{A vector with (starting) values for the response consistency parameters of each subject.}
\item{tol}{tolerance for the EM algorithm; estimation stops when the change in logLikelihood is smaller than tol.}
\item{maxiter}{maximum number of iterations of the EM algorithm.}
\item{A.est}{A list with each element a matrix with integers to fix (value 0) that parameter in A or freely estimate (value > 1) it. Giving different elements the same integer > 0 will contrain these elements to have the same value.}
\item{prior.est}{A vector with integers to fix or estimate the prior strategy probabilities. Numbering is similar to that for A.est.}
\item{b.est}{A vector with integers to fix or estimate the response consistency parameter. Numbering is similar to that for A.est}
\item{A.group}{a vector indicating for each subject which transition matrix (element in list A) applies.}
\item{verbose}{(logical) to indicate whether details of the EM iterations should be printed to screen.}
\item{b.min}{lower bound for the consistency parameter.}
}
\details{This function estimates the Strategy Switch model (Speekenbrink et al., in press). The Strategy Switch model is formulated as a hidden Markov model
in which the states are specific strategies of responding to cue patterns. Design matrices for the cue patterns are given in the list X. Each strategy specifies a predicted probability of responses for each cue pattern. These predictions are on a logit scale. The strategies are given as rows in the matrix Z, which contain coefficients such that X[i]][j,]%*%t(Z[k,]) results in the required logit value for subject i on trial j for strategy k. The example specifies the Z matrix for the Constant error strategies defined in Speekenbrink et al (2009). The are the Random, Singleton, Single Cue 1, Single Cue 2, Single Cue 3, Single Cue 4 and Multi-Cue strategies.
}
\value{An list with estimated A, b, prior, etc.}
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
\dontrun{mod <- StrategySwitch(y=y,X=X,Z=Z,prior=prior,A=A,b=2.8,prior.est=rep(0,7),A.est=TRUE,b.est=TRUE,tol=1e-5,maxiter=2000,A.group=c(1,2),verbose=T,b.min=0.4054651)
## compute BIC
-2*mod$LL + 2*mod$df}
}
