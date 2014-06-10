

# evol.model class

JiveMod <- list(lik.fun = function(),
			   hprior = function(),
			   win.size = list(),
			   init.val = list())

hpriorBM <- function(pars, list(dist) ){
	

}
	

# phy - simmap object (OUM) or phylo object (for BM or OU1)
makeJive <- function(phy, traits, model1="BM", model2="OU"){


	

	if (model == BM){
	  JiveMod <- list(lik.fun = likBM, # call to a function
			   hprior.fun = hpriorBM, # 
			   win.size = list(),
			   init.val = list())
	}
	iflse 

}


user <- list(id = 1,
             password = '41bfe136a536b7749104415bc364df2e',
             email = 'jmw@johnmyleswhite.com')
 
class(user) <- 'user'
 
user


jiveModel <- function(lik, winsize, initpar, hprior) {
	ll <- list(lik=function(){}, winsize=function(){}, initpar=function(){}, hprior=function(){})
	class(ll) <- "jiveModel"
	invisible(ll)
}

my.jiveModel <- jiveModel(likBM, initWinSize, initUpdateFreq, expMat)