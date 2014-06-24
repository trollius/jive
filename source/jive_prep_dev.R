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
make.jive <- function(simmap, traits, model_var="OU1", model_mean="BM", model_lik="Multinorm"){

	jive <- list()

	if (name.check(simmap, traits) != "OK") {
	
		print ("Species do not match")
	
	} else {

		if (model_var %in% c("BM", "OU1")) { # this is needed to overcome phytools limitation about making simmap obj from a trait with a single regime# td$phy$node.label <- rep("1", n - 1)
			jive$data$map <- matrix(rep(1, ((dim(traits)[1]) * 2) - 2))
		} else {
			jive$data$map <- relSim(simmap)$mapped.edge
		}
				
		jive$data$traits 					<- traits
		jive$data$counts 					<- apply(traits, 1, function (x) {sum( !is.na(x) )})
		jive$data$tree   					<- simmap
		jive$data$vcv    					<- vcv(simmap)
		jive$data$nreg   					<- dim(jive$data$map)[2]
		
		print(jive$data$nreg)
		
		
		if (model_lik == "Multinorm") {
			jive$lik$model 					<- likMultinorm
			jive$lik$mspws 					<- initWinSizeMVN(jive$data$traits)$msp
			jive$lik$sspws 					<- initWinSizeMVN(jive$data$traits)$ssp
			jive$lik$mspinit				<- initParamMVN(jive$data$traits)$mspA
			jive$lik$sspinit				<- initParamMVN(jive$data$traits)$sspA
				
		}
		
		if (model_mean == "BM" ) {
			jive$prior_mean$model 			<- likBM
			jive$prior_mean$init  			<- initParamMBM(jive$data$traits)  # check
			jive$prior_mean$ws	  			<- initWinSizeMBM(jive$data$traits)  # check
			jive$prior_mean$hprior$r		<- make.hpfun("Gamma", c(1.1,5))		  # sigma
			jive$prior_mean$hprior$m		<- make.hpfun("Uniform", c(-10000,10000)) # anc.mean

		
		}
		
		if (model_var == "BM" ) {
			jive$prior_var$model 			<- likBM
			jive$prior_var$init  			<- initParamVBM(jive$data$traits) # check
			jive$prior_var$ws	  			<- initWinSizeVBM(jive$data$traits) # check
			jive$prior_var$hprior$r			<- make.hpfun("Gamma", c(1.1,5)) # sigma
			jive$prior_var$hprior$m			<- make.hpfun("Gamma", c(1.1,5)) # anc.mean

		
		}

		
		if (model_var == "OU1" || model_var == "OUM") {
			jive$prior_var$model 			<- likOU
			jive$prior_var$init  			<- initParamVOU(jive$data$traits, jive$data$nreg)  # check
			jive$prior_var$ws	  			<- initWinSizeVOU(jive$data$traits, jive$data$nreg)  # check
			jive$prior_var$hprior$a			<- make.hpfun("Gamma", c(1.1,5)) # alpha
			jive$prior_var$hprior$r			<- make.hpfun("Gamma", c(1.1,5)) # sigma
			jive$prior_var$hprior$m			<- make.hpfun("Gamma", c(1.1,5)) # anc.mean
			for (i in 1:jive$data$nreg){									 # theta
				ti = paste("t",i,sep="")
				jive$prior_var$hprior[[ti]]			<- make.hpfun("Gamma", c(1.1,5))
			}
		
		}
		
	}	
	return(jive)

}




