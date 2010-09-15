#smc_resample_systematic <- function(logWeights,...) {
#    sw <- sum(weights)
#    u1 <- runif(1,0,sw/N)
#    u <- c(0,seq(sw/N,sw*(N-1)/N,length=N-1)) + u1
#    times <- hist(u,breaks=c(0,cumsum(weights)),plot=FALSE)$counts
#    ids <- rep(1:N,times=times)
#    return(ids)
#}

smc_resample_multinomial <- function(logWeights,...) {
    N <- length(logWeights)
    weights <- exp(logWeights)
	times <- rowSums(rmultinom(N,size=1,prob=weights))
    ids <- rep(1:N,times=times)
    return(ids)
}

smc_resample_residual <- function(logWeights,...) {
    N <- length(logWeights)
	weights <- exp(logWeights)
	Nw <- floor(N*weights/sum(weights))
	p <- weights - Nw/N
	times <- Nw + rowSums(rmultinom(N,size=1,prob=p))
    ids <- rep(1:N,times=times)
    return(ids)
}

setMethod("Resample",signature(object="ParticleMatrix"),
	function(object,type=c("systematic","residual","multinomial","stratified"),...) {
		type <- match.arg(type)
		logWeights <- getLogWeights(object,...)
		N <- object@N
		switch(type,
			systematic = {
				ids <- .Call("smc_resample_systematic",logWeights = logWeights)
			},
			residual = {
				ids <- smc_resample_residual(logWeights,...)
			},
			multinomial = {
				ids <- smc_resample_multinomial(logWeights,...)
			},
			stratified = {
				ids <- .Call("smc_resample_stratified",logWeights = logWeights)
			}
		)
		if(object@p_margin == 2) object@particles <- object@particles[,ids] else object@particles <- object@particles[ids,]
		object@logWeights <- rep(log(1/N),N)
		object@unifWeights <- TRUE
		object
	}
)

setMethod("doResample",signature(object="ParticleMatrix"),
	function(object,type=c("systematic","residual","multinomial","stratified"),...) {
        name <- deparse(substitute(object))
		type <- match.arg(type)
		weights <- getWeights(object,...)
		N <- object@N
				switch(type,
			systematic = {
				ids <- .Call("smc_resample_systematic",logWeights = logWeights)
			},
			residual = {
				ids <- smc_resample_residual(logWeights,...)
			},
			multinomial = {
				ids <- smc_resample_multinomial(logWeights,...)
			},
			stratified = {
				ids <- .Call("smc_resample_stratified",logWeights = logWeights)
			}
		)
		if(object@p_margin == 2) object@particles <- object@particles[,ids] else object@particles <- object@particles[ids,]
		object@logWeights <- rep(log(1/N),N)
		object@unifWeights <- TRUE
		assign(name,object,envir=parent.frame())
        return(invisible())
	}
)
