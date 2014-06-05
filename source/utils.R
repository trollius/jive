
# define classes for thermodynamic integration
DefineClasses <- function(ncategories, beta.param){ # if reviewers ask - see MrBayes, BEAST and Xia et al 2011 Sys Bio
     K <- ncategories-1
     k <- 0:K
     b <- k/K
     temperature <- rev(b^(1/beta.param))
}

# initialize update frequencies         
InitUpdateFreq <- function(jive.obj, update.freq=NULL){

	# init chain
	# calculate update frequencies
	if (!is.null(update.freq)) {
		update.freq	<- cumsum(update.freq/sum(update.freq))
	} else {
		update.freq	<- c(0.4,0.1,0.5)	
	}
	
	return(update.freq)

}


# user-defined window.sizes are not there - needs implementaiton
InitWinSize <- function(jive, window.sizes=NULL){
		
	# init window sizes
	ll.ws <- list()
	xx    <- apply(jive$traits, 1, sd, na.rm = T) # CHECK IF ITS SD OR VAR
	yy    <- sd(xx)
	
	ll.ws[["msp"]] <- xx 
	ll.ws[["ssp"]] <- xx
	ll.ws[["mvn"]] <- c(2*yy, yy) # anc.state of means windows size, evol rate of means window size
	
	if (jive$model == "BM") {
		ll.ws[["svn"]] <- c(2*yy, yy) # anc.state of sigmas windows size, evol rate of sgimas window size
	} else { 
		ll.ws[["svn"]] <- c(1.5,yy, 2*yy, rep(2*yy, jive$nreg)) # alpha from max likelihood on observed std dev	!!!!!################# IMPORTANT CNAnGE FOR ALPHA
	}
	
	return(ll.ws)
	
}
