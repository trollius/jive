
# PRIORS (they return log prior probabilities)
priorGamma <- function(i, a, b) {
	# Un-normalized gamma prior with shape parameters a, b.
	# 
	# Args:
	# 	i:  vector of parameters or a parameter
	#	a:  a parameter of a gamma distribution
	#	b:  b parameter of a gamma distribution
	# 
	# Returns:
	#	Vector of prior(s).

    (a - 1) * log(i) + (-b * i)
}

priorExponential <- function(i, l) {
	# Un-normalized exponential prior with mean l.
	# 
	# Args:
	# 	i:  vector of parameters or a parameter
	#	l:  shape parameter of an exponential distribution
	# 
	# Returns:
	#	Vector of prior(s).

     -(1./l) * i
}

priorUniform <- function(i, m, M) {
	# Uniform prior with boundaries m, M.
	# 
	# Args:
	#  	i:  vector of parameters or a parameter
	# 	m:  minimum boundary
	#	M:  maximum boundary
	# 
	# Returns:
	#	Vector of prior(s).
    # 
    if (i < m || i > M) {
         p <- -Inf
	} else {
         p <- -log(M - m)
    }
}
