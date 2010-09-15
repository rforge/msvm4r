\name{ParticleMatrix-class}

\docType{class}

\alias{ParticleMatrix-class}

\alias{getParticles,ParticleMatrix-method}
\alias{setParticles,ParticleMatrix-method}
\alias{ParticleMove,ParticleMatrix-method}
\alias{doParticleMove,ParticleMatrix-method}
\alias{SmcIterate,ParticleMatrix-method}
\alias{doSmcIterate,ParticleMatrix-method}
\alias{UpdateWeights,ParticleMatrix-method}
\alias{doUpdateWeights,ParticleMatrix-method}
\alias{ESS,ParticleMatrix-method}
\alias{getLogWeights,ParticleMatrix-method}
\alias{getWeights,ParticleMatrix-method}
\alias{getNormWeights,ParticleMatrix-method}

\title{Class "ParticleMatrix"}

\description{Basic class to implement a Sequential Monte Carlo sampler.}

\section{Slots}{

	\describe{
    \item{\code{Particles}:}{Matrix with particle values}
    
    \item{code{p_margin}:}{Integer indicating which margin in the matrix represents particle ids}

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
    \item{getParticles}{\code{signature(object = "ParticleMatrix")}: return particles}
    \item{setParticles}{\code{signature(object = "ParticleMatrix")}: set particles}
    \item{ParticleMove}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{doParticleMove}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{SmcIterate}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{doSmcIterate}{\code{signature(object = "ParticleMatrix")}: move particles}
    
    \item{UpdateWeights}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{doUpdateWeights}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{ESS}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{getLogWeights}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{setLogWeights}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{getWeights}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{getNormWeights}{\code{signature(object = "ParticleMatrix")}: move particles}

  }

}

\author{Maarten Speekenbrink}

\keyword{classes}
