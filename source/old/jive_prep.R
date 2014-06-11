library(geiger)

## tranform simmap into map, input - simmap object
relSim <- function(x) {
	
	foo<-function(x) {
		x/sum(x)
	}
	
	x$mapped.edge <- t(apply(x$mapped.edge, 1, FUN=foo))
	x$mapped.edge <- x$mapped.edge[, order(colnames(x$mapped.edge))]
	return(x)
				
}

# traits - matrix, where lines are vector of observations for each species, with NA for no data.
# phy - simmap object (OUM) or phylo object (for BM or OU1)
makeJive <- function(phy, traits, model="BM"){

	if (class(phy) != "phylo") {
		if (model %in% c("BM", "OU1")) {
			stop("'phy' should be an object of class 'phylo'")
		} else {
			stop("'phy' should be a 'simmap' object")
		}
	}
	
	if (class(traits) != "matrix") {
		stop("'traits' should be an object of class 'matrix'")
	}
	
	if (!(model %in% c("BM", "OU1", "OUM"))) {
		stop("Only BM, OU1, and OUM models are supported for now.")
	}
	
	ll <- list()
	td <- treedata(phy, traits)
	
	if (model %in% c("BM", "OU1")) { # this is needed to overcome phytools limitation about making simmap obj from a trait with a single regime
		# td$phy$node.label <- rep("1", n - 1)
		ll[['map']] <- matrix(rep(1, (n * 2) - 2))
	} else {
		ll[['map']] <- relSim(td$phy)
	}
			
	ll[['traits']] <- td$data
	ll[['counts']] <- apply(td$data, 1, function (x) {sum( !is.na(x) )})
	ll[['tree']]   <- td$phy
	ll[['vcv']]    <- vcv(td$phy)
	ll[['model']]  <- model
	ll[['nreg']]   <- dim(td$phy)[2]
	
	
	return(ll)
	
}


constrBM <- function(data){
	
	ll <- list()
	pp <- initWinSize(data)
	dd <- initParam(data)
	
	ll$lik	<- likBM # function call
	ll$hp	<- hyperpBM # function call
	ll$ws	<- pp # list
	ll$init	<- dd # list


	return(ll)

}


constrOU1 <- function(){



}


constrOUM <- function(){



}
