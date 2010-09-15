setMethod("ESS",signature(object="ParticleBase"),
	function(object,...) {
		#w <- exp(object@logWeights)
		sumw <- sum(exp(object@logWeights))
		sumsq <- sum(exp(2*object@logWeights))
		exp(-log(sumsq) + 2*log(sumw))
	}
)
