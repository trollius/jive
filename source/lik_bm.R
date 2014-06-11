
#conditional prior on sd of species-specific normal likelihoods under GBM model, x - vector of values
likBM<-function(pars, x, tree){#M - ancestral mean, S - trend, S0 - starting point of a trend, ti - total phylogenetic time, sig.sq  - sigma^2 (phylogenetic variance)
	
	vcv.m  <- vcv(tree)
	Y      <- as.matrix(x)
	n      <- dim(vcv.m)[1]
	m      <- matrix(1, n, 1)
	m[, ]  <- pars[3] # ancestral mean
	sig.sq <- pars[2] # sigma
	DET    <- determinant(sig.sq * vcv.m, logarithm=T)

	log.lik.GBM <- try((-n/2 * log(2 * pi) - (as.numeric(DET$modulus))/2 - 1/2 * (t(Y - m)%*%ginv(sig.sq * vcv.m)%*%(Y - m))), silent=T)
	
	if (is.na(log.lik.GBM) | (class(log.lik.GBM) == "try-error" )) {
			return(-Inf)
	} else {
		return(log.lik.GBM)
	}

}

	