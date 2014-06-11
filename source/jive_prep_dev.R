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




