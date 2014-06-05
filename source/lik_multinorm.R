# input: means = vector of means, sigmas = vector of sigmas,SP_DATA = matrix of species observations, N_OBSERV vector of observation counts
# does: calculate individual log-likelihoods for each species based on normal distribution
LikMultinorm <- function(means, sigmas, traits, counts){#m - mean (horizontal), s - sigma^2 (horizontal), vec - observations for a species
	
	
	log.lik.MN <- -counts/2 * log(2 * pi) - 1/2 * counts * log(sigmas) - 1/2 * (apply((traits - means)^2, 1, sum, na.rm=T)/sigmas)
	
	if (is.na(sum(log.lik.MN))) {
			return(-Inf)
	} else {
		return(log.lik.MN)
	}

	
	
	
} # Gaussian density

