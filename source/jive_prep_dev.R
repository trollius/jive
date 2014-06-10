library(geiger)
# jive structure
# jive$data			# $traits 
					# $counts 
					# $tree  
					# $vcv 
					# $map
					# $nreg					
				
# jive$lik			# $model
					# $wsmsp
					# $wsssp
					# $initmsp
					# $initssp
				
# jive$prior_mean 	# $model
					# $init
					# $ws
					# $hprior
						
# jive$prior_var	# $model
					# $init
					# $ws
					# $hprior

					
# traits - matrix, where lines are vector of observations for each species, with NA for no data.
# phy - simmap object (OUM) or phylo object (for BM or OU1)
make.jive <- function(simmap, traits, model_var="OU1", model_mean="BM", model_lik="MVN"){

	jive <- list()

	if (name.check(simmap, traits) != "OK") {
	
		print ("Species do not match")
	
	} else {

		if (model_var %in% c("BM", "OU1")) { # this is needed to overcome phytools limitation about making simmap obj from a trait with a single regime# td$phy$node.label <- rep("1", n - 1)
			jive$data$map <- matrix(rep(1, ((dim(traits)[1]) * 2) - 2))
		} else {
			jive$data$map <- relSim(simmap)
		}
				
		jive$data$traits 					<- traits
		jive$data$counts 					<- apply(traits, 1, function (x) {sum( !is.na(x) )})
		jive$data$tree   					<- simmap
		jive$data$vcv    					<- vcv(simmap)
		jive$data$nreg   					<- dim(jive$data$map)[2]
		
		print(jive$data$nreg)
		
		
		if (model_lik == "MVN") {
			jive$lik$model 					<- likMVN
			jive$lik$mspws 					<- initWinSizeMVN(jive$data$traits)$msp
			jive$lik$sspws 					<- initWinSizeMVN(jive$data$traits)$ssp
			jive$lik$mspinit				<- initParamMVN(jive$data$traits)$mspA
			jive$lik$sspinit				<- initParamMVN(jive$data$traits)$sspA
				
		}
		
		if (model_mean == "BM") {
			jive$prior_mean$model 			<- likBM
			jive$prior_mean$init  			<- initParamBM(jive$data$traits)
			jive$prior_mean$ws	  			<- initWinSizeBM(jive$data$traits)
			jive$prior_mean$hprior_m		<- make.hpfun("Uniform", c(-10000,10000))
			jive$prior_mean$hprior_r		<- make.hpfun("Gamma", c(1.1,5))
		
		}

		
		if (model_var == "OU1" || model_var == "OUM") {
			jive$prior_var$model 			<- likOU
			jive$prior_var$init  			<- initParamOU(jive$data$traits, jive$data$nreg)
			jive$prior_var$ws	  			<- initWinSizeOU(jive$data$traits, jive$data$nreg)
			jive$prior_var$hprior_a			<- make.hpfun("Gamma", c(1.1,5))
			jive$prior_var$hprior_m			<- make.hpfun("Gamma", c(1.1,5))
			jive$prior_var$hprior_r			<- make.hpfun("Gamma", c(1.1,5))
			jive$prior_var$hprior_t			<- make.hpfun("Gamma", c(1.1,5))
		
		}
	}	
	return(jive)

}

##-------------------------- set hyperpriors
make.hpfun <-function(hpf="Uniform", hp.pars){


		if (hpf == "Uniform"){
			my.f <- function(x){
				hp <- sum(dunif(x, min=hp.pars[1], max=hp.pars[2], log=TRUE))
				return(hp)
			}
			
		}
		
		if (hpf == "Gamma"){
			my.f <- function(x){
				hp <- sum(dgamma(x, shape=hp.pars[1], scale=hp.pars[2], log=TRUE))
				return(hp)
			}
		}
		
		if (hpf == "Normal"){
			my.f <- function(x){
				hp <- sum(dnorm(x, mean=hp.pars[1], sd=hp.pars[2], log=TRUE))
				return(hp)
			}
		}	
		
		
	return(my.f)
}

## tranform simmap into map, input - simmap object
relSim <- function(x) {
	
	foo<-function(x) {
		x/sum(x)
	}
	
	x$mapped.edge <- t(apply(x$mapped.edge, 1, FUN=foo))
	x$mapped.edge <- x$mapped.edge[, order(colnames(x$mapped.edge))]
	return(x)
				
}

##-------------------------- initialize windows sizes functions
initWinSizeMVN <- function (x){

	ws		<- list()
	xx		<- apply(x, 1, sd, na.rm = T) # CHECK IF ITS SD OR VAR
	ws$msp	<- xx 
	ws$ssp	<- xx
	return(ws)

}
 
# input is trait matrix, rows are species, cols - observations
initWinSizeBM <- function(x){
	
	xx <- sd(apply(x, 1, sd, na.rm = T)) # CHECK IF ITS SD OR VAR
	ws <- c(2 * xx, xx) # anc.state of sigmas windows size, evol rate of sgimas window size

	return(ws)
	
}

# input is trait matrix, rows are species, cols - observations
initWinSizeOU <- function(x, nreg){
	
	xx <- sd(apply(x, 1, sd, na.rm = T)) # CHECK IF ITS SD OR VAR
	ws <- c(1.5, xx, 2 * xx, rep(2 * xx, nreg)) # alpha from max likelihood on observed std dev <------------------------ alpha parameter to adjust
	return(ws)
	
}

##-------------------------- initialize start parameters values functions
# initialize MCMC parameters
initParamMVN <- function (x){

	init  <- list()
	init$mspA  <- apply(x, 1, mean, na.rm = T) # initialize means for species
	init$sspA  <- apply(x, 1, var, na.rm = T) # initialize sigma.sq for species
	
	return(init)
}

# initialize MCMC parameters				   
initParamBM <- function(x){
		
	# initialize MCMC parameters
	init <- c(mean(apply(x, 1, mean, na.rm = T)), var(apply(x, 1, mean, na.rm = T))) # could be a random number, initialize mean for MVN
	return(init)

}				   

# initialize MCMC parameters
initParamOU <- function(x, nreg){
		
	init <- c(runif((nreg+3), 0.5, 3)) # could be aither more realistic values such as means and sds of true data
	#init <- c(2.941516,2.139533,1.299683,1.364224) just a check
	return(init)

}	

