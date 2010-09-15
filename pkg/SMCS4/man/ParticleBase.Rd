\name{ParticleBase-class}

\docType{class}

\alias{ParticleBase-class}

\alias{getParticles,ParticleBase-method}
\alias{setParticles,ParticleBase-method}
\alias{ParticleMove,ParticleBase-method}
\alias{doParticleMove,ParticleBase-method}
\alias{SmcIterate,ParticleBase-method}
\alias{doSmcIterate,ParticleBase-method}
\alias{UpdateWeights,ParticleBase-method}
\alias{doUpdateWeights,ParticleBase-method}
\alias{ESS,ParticleBase-method}
\alias{getLogWeights,ParticleBase-method}
\alias{getWeights,ParticleBase-method}
\alias{getNormWeights,ParticleBase-method}

\title{Class "ParticleBase"}

\description{Basic class to implement a Sequential Monte Carlo sampler.}

\section{Slots}{

	\describe{

		logWeights = "vector", # log of particle weights
		unifWeights = "logical", # are current weights uniform?
		p_move = "function",
		lW_update = "function",
		resampleC = "numeric",
		N = "integer"

	\item{\code{logWeights}:}{Vector containing the log (unnormalised) particle weights.}

	\item{\code{unifWeights}}{Logical indicating whether the logWeights are uniform.}

	\item{\code{p_move}:}{Function to move the particles to a new position.}
	
	\item{\code{mcmc_move}:}{Function to perform a Monte Carlo Markov Chain move.}

	\item{\code{lW_update}:}{Function to do update the \code{logWeights}.}

	\item{\code{resampleC}:}{Numeric value (between 0 and 1) indicating when to perform resampling.}

	}
}

\section{Methods}{
\describe{
    \item{getParticles}{\code{signature(object = "ParticleBase")}: return particles}
    \item{setParticles}{\code{signature(object = "ParticleBase")}: set particles}
    \item{ParticleMove}{\code{signature(object = "ParticleBase")}: move particles}
    \item{doParticleMove}{\code{signature(object = "ParticleBase")}: move particles}
    \item{SmcIterate}{\code{signature(object = "ParticleBase")}: move particles}
    \item{doSmcIterate}{\code{signature(object = "ParticleBase")}: move particles}
    
    \item{UpdateWeights}{\code{signature(object = "ParticleBase")}: move particles}
    \item{doUpdateWeights}{\code{signature(object = "ParticleBase")}: move particles}
    \item{ESS}{\code{signature(object = "ParticleBase")}: move particles}
    \item{getLogWeights}{\code{signature(object = "ParticleBase")}: move particles}
    \item{setLogWeights}{\code{signature(object = "ParticleBase")}: move particles}
    \item{getWeights}{\code{signature(object = "ParticleBase")}: move particles}
    \item{getNormWeights}{\code{signature(object = "ParticleBase")}: move particles}

  }

}

\author{Maarten Speekenbrink}

\keyword{classes}
