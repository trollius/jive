library(geiger)

## tranform simmap into map, input - simmap object
RelSim <- function(x) {
	
	foo<-function(x) {
		x/sum(x)
	}
	
	x$mapped.edge <- t(apply(x$mapped.edge, 1, FUN=foo))
	x$mapped.edge <- x$mapped.edge[, order(colnames(x$mapped.edge))]
	return(x)
				
}

# traits - matrix, where lines are vector of observations for each species, with NA for no data.
# phy - simmap object (OUM) or phylo object (for BM or OU1)
MakeJive <- function(phy, traits, model="BM"){
	
	ll <- list()
	td <- treedata(phy, traits)
	
	if (model %in% c("BM", "OU1")) { # this is needed to overcome phytools limitation about making simmap obj from a trait with a single regime
		# td$phy$node.label <- rep("1", n - 1)
		ll[["map"]]    <- matrix(rep(1, (n * 2) - 2))
	} else {
		ll[["map"]]    <- RelSim(td$phy)
	}
			
	ll[["traits"]] <- td$data
	ll[["counts"]] <- apply(td$data, 1, function (x) {sum( !is.na(x) )})
	ll[["tree"]]   <- td$phy
	ll[["vcv"]]    <- vcv(td$phy)
	ll[["model"]]  <- model
	ll[["nreg"]]   <- dim(td$phy)[2]
	
	return(ll)
	
}




