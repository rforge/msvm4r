setClass("ParticleMatrix",
	contains="ParticleBase",
	representation(
		particles = "matrix",
        p_margin = "integer"
	),
	prototype=prototype(
        p_margin = as.integer(1)
    )

)

setMethod("ParticleMean",signature(object="ParticleMatrix"),
	function(object,...) {
	  if(!object@unifWeights) {
	    w <- getNormWeights(object,...)
	    if(object@p_margin == 2) return(colSums(w*t(object@particles))) else return(colSums(w*object@particles))
	  } else {
	    if(object@p_margin == 2) return(rowMeans(object@particles)) else return(colMeans(object@particles))
	  }
	}
)

setMethod("ParticleCov",signature(object="ParticleMatrix"),
	function(object,...) {
	  if(!object@unifWeights) {
	    w <- getWeights(object,...)
	    if(object@p_margin == 2) return(cov.wt(t(object@particles),...)$cov) else return(cov.wt(object@particles,...)$cov)
	  } else {
	    if(object@p_margin == 2) return(cov(t(object@particles),...)) else return(cov(object@particles,...))
	  }
	}
)


setMethod("ParticleVar",signature(object="ParticleMatrix"),
	function(object,...) {
	  if(!object@unifWeights) {
	    w <- getWeights(object,...)
	    if(object@p_margin == 2) return(diag(cov.wt(t(object@particles),...)$cov)) else return(diag(cov.wt(object@particles,...)$cov))
	  } else {
	    if(object@p_margin == 2) return(apply(object@particles,1,var,...)) else return(apply(object@particles,2,var,...))
	  }
	}
)

