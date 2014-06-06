

defClasses <- function(ncat=10, beta.param=0.3){ 
	# Defines classes for thermodynamic integration.
	# For details of the method see Xia et al 2011 Sys Bio.
	#
	# Args:
	# 	ncategories: number of classes that will be used in thermodynamic integration.
	#	beta.param:  parameter describing the shape of a beta distribution.
	# 
	# Returns:
	#	The vector of temperatures for thermodynamic integration.
	
	
    K <- ncat-1
    k <- 0:K
    b <- k/K
    temp<- rev(b^(1/beta.param))
	return(temp)
}

       
initUpdateFreq <- function(update.freq=NULL){
	# Initializes update frequencies for likelihood and two prior levels.
	#
	# Args:
	# 	update.freq: the vector (length = 3) of update frequencies (likelihood, priorMBM, priorVOU/VBM).
	#
	# Returns:
	#	The vector of update frequencies which sums to 1. 
	
	if (length(update.freq) != 3 && !is.null(update.freq)) {
		stop("Update.freq must contain 3 elements" )
	}
	
	
	# calculate update frequencies
	if (!is.null(update.freq)) {
		update.freq	<- cumsum(update.freq/sum(update.freq))
	} else {
		update.freq	<- c(0.4,0.1,0.5)	
	}
	
	return(update.freq)

}


initWinSize <- function(jive, window.sizes=NULL){
	# Calculates window sizes for sampling proposals.
	# User-defined window sizes are not supported at this stage.
	#
	# Args:
	# 	jive: jive.object (see makeJive function)
	#
	# Returns:
	#   A list of windows sizes for: msp - species-specific means, 
	#							     ssp - species-specific variances,
	#								 mvn - prior on means (MBM),
	#								 svn - prior on variances (VBM/VOU).
	
	ll.ws <- list()
	xx    <- apply(jive$traits, 1, sd, na.rm = T) # CHECK IF ITS SD OR VAR
	yy    <- sd(xx)
	
	ll.ws[["msp"]] <- xx 
	ll.ws[["ssp"]] <- xx
	ll.ws[["mvn"]] <- c(2 * yy, yy) # anc.state of means windows size, evol rate of means window size
	
	if (jive$model == "BM") {
		ll.ws[["svn"]] <- c(2 * yy, yy) # anc.state of sigmas windows size, evol rate of sgimas window size
	} else { 
		ll.ws[["svn"]] <- c(1.5, yy, 2 * yy, rep(2 * yy, jive$nreg)) # alpha from max likelihood on observed std dev <------------------------ alpha parameter to adjust
	}
	
	return(ll.ws)
	
}
