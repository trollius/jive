

# UPDATE PARAMETERS
# sliding window
SlidingWin <- function(i, d, M) {
	# i = current value; d = window size; M = max allowed value 
	ii <- abs(i + (runif(length(i), 0, 1) - 0.5) * d) #MrBayes trick
	# reflection at 0 and M
	#if (ii>M) {ii=abs((M-(ii-M)))}
	#if (ii>M) {ii=i}
	return(ii)
}

# allowing negative values for proposal	
SlidingWinUnconst <- function(i, d) {
     # i = current value; d = window size; M = max allowed value 
     ii <- i + (runif(length(i), 0, 1) - 0.5) * d #MrBayes trick
     ii
}     
     
# ronquist multiplier
MultiplierProposal <- function(i, d, u) {
     lambda <- 2 * log(d)
     m <- exp(lambda * (u - 0.5))
     ii <- i * m
     ii
}

# PRIORS (they return log prior probabilities)
PriorGamma <- function(i, a, b) {
     # un-normalized gamma with shape parameters a, b
     # i parameter (or vector of parameters)
     (a - 1) * log(i) + (-b * i)
}

PriorExponential <- function(i, l) {
     # un-normalized exponential with mean l
     # i parameter (or vector of parameters)
     -(1./l) * i
}

PriorUniform <- function(i, m, M) {
     # uniform prior with boundaries m, M
     if (i < m || i > M) {
          p <- -Inf
	 } else {
          p <- -log(M - m)
     }
}
