setClass("ParticleBase",
	representation(
		logWeights = "vector", # log of particle weights
		unifWeights = "logical", # are current weights uniform?
		p_move = "function",
		mcmc_move = "function",
		lW_update = "function",
		logLik = "function",
		resampleC = "numeric",
		N = "integer"
	)
)


setMethod("getParticles",signature(object="ParticleBase"),
    function(object,...) {
        return(object@particles)
    }
)

setReplaceMethod(
    f = "setParticles",
    signature = "ParticleBase",
    definition = function(object,value) {
        object@particles <- value
        return(object)
    }
)

setMethod("ParticleMove",signature(object="ParticleBase"),
	function(object,...) {
		object@particles <- object@p_move(getParticles(object,...),...)
		object
	}
)

setMethod("doParticleMove",signature(object="ParticleBase"),
	function(object,...) {
        name <- deparse(substitute(object))
        name <- paste(name,"@particles",sep="")
		#setParticles(object,...) <- object@p_move(object,...)
		assign(name,object@p_move(getParticles(object,...),...),envir=parent.frame())
        return(invisible())
	}
)

setMethod("McmcMove",signature(object="ParticleBase"),
	function(object,...) {
		object@particles <- object@mcmc_move(getParticles(object,...),...)
		object
	}
)

setMethod("doMcmcMove",signature(object="ParticleBase"),
	function(object,...) {
        name <- deparse(substitute(object))
        name <- paste(name,"@particles",sep="")
		#setParticles(object,...) <- object@p_move(object,...)
		assign(name,object@mcmc_move(getParticles(object,...),...),envir=parent.frame())
        return(invisible())
	}
)


setMethod("SmcIterate",signature(object="ParticleBase"),
	function(object,...) {
		object <- ParticleMove(object,...)
		object <- UpdateWeights(object,...)
		ess <- ESS(object,...)
		if(ess < (object@resampleC*object@N)) object <- Resample(object,...)
		object
	}
)

setMethod("doSmcIterate",signature(object="ParticleBase"),
	function(object,...) {
        name <- deparse(substitute(object))
		doParticleMove(object,...)
        doUpdateWeights(object,...)
        ess <- ESS(object,...)
		if(ess < (object@resampleC*object@N)) doResample(object,...)
		assign(name,object,envir=parent.frame())
        return(invisible())
	}
)

setMethod("UpdateWeights",signature(object="ParticleBase"),
	function(object,data,...) {
		logWeights <- object@lW_update(getParticles(object,...),getLogWeights(object,...),...)
		# normalize to sensible values
		max <- max(-.Machine$double.xmax,logWeights)
		object@logWeights <- logWeights - max
		
		#dMaxWeight = -std::numeric_limits<double>::infinity();
    # for(int i = 0; i < N; i++)
    #   dMaxWeight = max(dMaxWeight, pParticles[i].GetLogWeight());
    # for(int i = 0; i < N; i++)
    #   pParticles[i].SetLogWeight(pParticles[i].GetLogWeight() - (dMaxWeight));  
		
		object@unifWeights <- FALSE
		object
	}
)

setMethod("doUpdateWeights",signature(object="ParticleBase"),
	function(object,data,...) {
        name <- deparse(substitute(object))
        name <- paste(name,"@logWeights",sep="")
		logWeights <- object@lW_update(getParticles(object,...),getLogWeights(object,...),...)
		# normalize to sensible values
		max <- max(-.Machine$double.xmax,logWeights)
        logWeights <- logWeights - max
        assign(paste(name,"@logWeights",sep=""),logWeights,envir=parent.frame())
		#object@unifWeights <- FALSE
		assign(paste(name,"@unifWeights",sep=""),FALSE,envir=parent.frame())
        return(invisible())
	}
)

setMethod("ESS",signature(object="ParticleBase"),
	function(object,...) {
		#w <- exp(object@logWeights)
		sumw <- sum(getWeights(object,...))
		sumsq <- sum(exp(2*getLogWeights(object,...)))
		exp(-log(sumsq) + 2*log(sumw))
	}
)

setMethod("getLogWeights",signature(object="ParticleBase"),
  function(object,...) {
    return(object@logWeights) 
  }
)

setReplaceMethod(
    f = "setLogWeights",
    signature = "ParticleBase",
    definition = function(object,value) {
        object@logWeights <- value
        return(object)
    }
)

setMethod("getWeights",signature(object="ParticleBase"),
  function(object,...) {
    if(!object@unifWeights) return(exp(getLogWeights(object,...))) else return(rep(1,object@N))
  }
)

setMethod("getNormWeights",signature(object="ParticleBase"),
  function(object,...) {
    if(!object@unifWeights) {
        weights <- getWeights(object,...)
        weights <- weights/sum(weights)
        return(weights)
    } else {
        return(rep(1/object@N,object@N))
    }
  }
)

setMethod("getN",signature(object="ParticleBase"),
    function(object,...) {
        return(object@N)
    }
)

setMethod("logLik",signature(object="ParticleBase"),
    function(object,...) {
        return(object@logLik(getParticles(object,...),...))
    }
)

