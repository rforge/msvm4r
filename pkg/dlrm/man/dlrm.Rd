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
filter.only=FALSE,verbose=FALSE,criterion=c("logLike","parameter"),
method="BFGS",hessian=FALSE,switch.LL=.5,switch.wait=5)
}
\arguments{
\item{formula}{a formula specifying the regression model (coefficients vary over
time).}
\item{data}{(optional) data framem to evaluate the formula.}
\item{maxit}{maximum number of iterations of the EM algorithm.}
\item{ws}{a vector with (starting values) for the coefficients (defaults to 0)}
\item{Sigma}{a matrix specyfing (starting values) of the covariance matrix of
initial coefficients.}
\item{A}{a matrix with (starting) values for the transition matrix of regression
coefficients.}
\item{Q}{a matrix with (starting) values for the covariance matrix of changes in 
the regression coefficients.}
\item{R}{a matrix with (starting) values for the observational error covariance
matrix.}
\item{Q.c}{a matrix}
\item{Sigma.c}{a matrix}
\item{ntimes}{a numeric vector specifying the length of each repeated time 
series.}
\item{tol}{tolerance for the EM algorithm}
\item{est.ws}{(logical) estimate the values in ws?}
\item{est.Sigma}{(logical) estimate the values in Sigma?}
\item{est.A}{(logical) estimate the values in A?}
\item{est.Q}{(logical) estimate the values in Q?}
\item{est.R}{(logical) estimate the values in R?}
\item{filter.only}{(logical) if TRUE, call returns the filtering estimates and 
does not estimate the (hyper-)parameters.}
\item{verbose}{(logical) print information to the terminal about the interations
in the EM algorithm?}
\item{criterion}{stopping criterion for the EM algorithm.}
\item{method}{method used in numerical optimization.}
\item{hessian}{(logical) return Hessian?}
\item{switch.LL}{value of logLik difference used when switching from EM to
numerical optimization.}
\item{switch.wait}{number of iterations before re-allowing numerical 
optimization.}
}
\details{This function estimates a dynamic linear regression model.
}
\value{A list with estimated A, b, prior, etc.}
\author{Maarten Speekenbrink}
\references{
Speekenbrink, M. & Shanks, D. R. (2009) Learning in a changing environment.
}
\examples{
## open StockMarket Prediction Task data
data(SMPT)
}
