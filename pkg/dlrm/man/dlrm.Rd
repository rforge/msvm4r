\name{dlrm}
\alias{dlrm}
\title{Dynamic linear regression model}
\description{
Estimates a dynamic linear regression model by the Kalman filter/smoother and
EM algorithm for parameter estimation.
}
\usage{
dlrm(formula,data,maxit=100,ws,Sigma,A,Q,R,
Q.c=NULL,Sigma.c=Q.c,ntimes=NULL,tol=1e-5,
est.ws=TRUE,est.Sigma=TRUE,est.A=TRUE,est.Q=TRUE,est.R=TRUE,
filter.only=FALSE,verbose=FALSE,criterion=c("logLik","parameter"),
method="BFGS",hessian=FALSE,switch.LL=.5,switch.wait=5)
}
\arguments{
\item{formula}{a formula specifying the regression model (regression 
coefficients vary over time).}
\item{data}{(optional) data frame to evaluate the formula.}
\item{maxit}{maximum number of iterations of the EM algorithm.}
\item{ws}{a vector with (starting values) for the regression coefficients 
(defaults to 0)}
\item{Sigma}{a matrix specyfing (starting values) of the covariance matrix of
initial coefficients.}
\item{A}{a matrix with (starting) values for the transition matrix of regression
coefficients.}
\item{Q}{a matrix with (starting) values for the covariance matrix of changes in 
the regression coefficients.}
\item{R}{a matrix with (starting) values for the observational error covariance
matrix.}
\item{Q.c}{a binary matrix specifying whether elements of Q are freely estimable 
(value 1) or fixed (value 0)}
\item{Sigma.c}{a binary matrix specifying whether elements of Sigma are freely 
estimable (value 1) or fixed (value 0)}
\item{ntimes}{a numeric vector specifying the length of each repeated time 
series.}
\item{tol}{tolerance level for the EM algorithm,}
\item{est.ws}{(logical) estimate the values in ws?}
\item{est.Sigma}{(logical) estimate the values in Sigma?}
\item{est.A}{(logical) estimate the values in A?}
\item{est.Q}{(logical) estimate the values in Q?}
\item{est.R}{(logical) estimate the values in R?}
\item{filter.only}{(logical) if TRUE, call returns the filtering estimates and 
does not estimate the (hyper-)parameters.}
\item{verbose}{(logical) print information to the terminal about the iterations
in the EM algorithm?}
\item{criterion}{stopping criterion for the EM algorithm. If argument is 
"logLik" the algorithm stops if the abolute change in the log likelihood is 
smaller than argument "tol". If the argument is "parameter", the algorithm stops
when the absolute change is smaller than argument "tol" for all estimated 
parameters.}
\item{method}{method used in numerical optimization. See "optim" for details.}
\item{hessian}{(logical) return the Hessian matrix of the parameter estimates?}
\item{switch.LL}{value of logLik difference used when switching from EM to
numerical optimization.}
\item{switch.wait}{number of iterations before re-allowing numerical 
optimization.}
}
\details{This function estimates a dynamic linear regression model via a 
combination of Expectation-Maximization (EM) and numerical optimization of the 
likelihood.

The model is for the (scalar) observations is

\eqn{y_t = \vec{x}_t \vec{w}_t + e_t}{y[t] = x[t] * w[t] + e[t]}

where \eqn{e_t}{e[t]} is a Normally distributed variable with mean 0 and 
variance R. The model for the regression coefficients is

\eqn{\vec{w}_{t} = \vec{w}_{t-1} + \vec{d}_t}{w[t] = w[t-1] + d[t]}

where \eqn{d_t}{d[t]} is a (multivariate) Normally distributed variable with
mean 0 and covariance matrix Q.

The function computes the posterior distribution p(w_{1:T}|y_{1:T},x_{1:T}) by 
the Kalman smoother. 
}

\value{A list with estimated A, b, prior, etc.}
\author{Maarten Speekenbrink}
\references{
Speekenbrink, M. & Shanks, D. R. (2010) Learning in a changing environment. 
Journal of Experimental Psychology: General.
}
\examples{
## open StockMarket Prediction Task data
data(SMPT)
mod <- dlrm(r ~ c1 + c2 - 1,data=SMPT,ws=c(0,0),Sigma=10*diag(2),
est.Sigma=FALSE,est.ws=FALSE,est.A=FALSE,Q=.005*diag(2),Q.c=diag(2),
ntimes=c(300,300),verbose=TRUE)
plot(mod$weight[1:300,1],type="l")
lines(mod$weight[301:600,1],col="red")

}
