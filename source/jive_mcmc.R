### MCMC TEMPLATE ###

# LOAD PACKAGES
library(OUwie)
library(nloptr)
library(ape)
library(MASS)
library(phytools)

source("lik_mvn.R")
source("lik_bm.R")
source("lik_ou.R")
source("lik_multinorm.R")
source("hpriors.R")
source("proposals.R")
source("utils.R")
source("jive_prep.R")


# MCMC MAIN ALGORITHM
jiveMCMC <- function(jive, log.file="jive_mcmc.log", sampling.freq=1000, print.freq=100, 
				ncat=1, beta.param=0.3, ngen=5500000, burnin=0, update.freq=NULL)
{

	# get data
		traits <- jive$traits
		model  <- jive$model
		counts <- jive$counts
		tree   <- jive$tree
		map    <- jive$map
				 
	# counter for the total chain length 
	real.iter <- 1 

	# here we define the chain length for each category
	it <- (ngen - burnin)/ncat 
		 if (ncat > 1) {
			  # here we get the heating parameter for a chain - scaling classes (MIGRATE)
			  beta.class <- defClasses(ncat, beta.param) 
		 } else {
			  beta.class <- 1
		 }
		 # start looping over scaling classes 
		 for (i.beta in 1:length(beta.class)){ 
			  # apply burning only to a first class
			  if (i.beta > 1) { 
				   burnin <- 0
			  }
			  #exctract temperature
			  temperature <- beta.class[i.beta]
				   acc <- 0 #acceptance count
				   for (iteration in 1:(it + burnin)) {
						hasting.ratio <- 0
						if (real.iter == 1){
							 # initialize update frequencies
							 update.freq <- initUpdateFreq(update.freq)
							 
							 # initialize window sizes
							 ws <- initWinSize(jive)
													 
							 # initialize MCMC parameters
							 mspA  <- apply(traits, 1, mean, na.rm = T) # initialize means for species
							 sspA  <- apply(traits, 1, var, na.rm = T) # initialize sigma.sq for species
							 mmvnA <- mean(mspA) # could be a random number, initialize mean for MVN
							 smvnA <- var(mspA) # could be a random number, initialize sigma.sq for MVN
							 bmouA <- c(runif(length(ws$svn), 0.5, 3)) # could be aither more realistic values such as means and sds of true data
							 
						}
						
						msp  <- mspA
						ssp  <- sspA
						mmvn <- mmvnA
						smvn <- smvnA
						bmou <- bmouA
						r    <- runif(1)
						
						# Update on all parameters level

						# update msp, ssp
						if (r < update.freq[1]) {
						
							 # 5 is just for now number of species updated
							 ind <- sample(1:length(msp), 5, replace=F)
							
							 if (runif(1) < 0.5) {
								  # updating random 5 values from the vector of means 
								  msp[ind] <- slidingWin(mspA[ind], ws$msp[ind]) 
							 } else {
								  # updating random 5 values from the vector of sigmas
								  ssp[ind] <- abs(slidingWin(sspA[ind], ws$ssp[ind]))
							 }
						}
						# update MVN parameters
						else if (r < update.freq[2]) {
						
							 if (runif(1) < 0.5){
								  # updating mmvn - negative values allowed
								  mmvn <- slidingWin(mmvnA, ws$mvn[1]) 
							 } else {
								  # updating smvn
								  smvn <- abs(slidingWin(smvnA, ws$mvn[2]) )
							 }
							 
						} else {# update BM parameters 
												 
							 ind <- sample(1:length(ws$svn), 1)
							 bmou[ind] <- abs(slidingWin(bmouA[ind], ws$svn[ind])) # updating bmou parameters
							 
						}
						

						# Hyperpriors on all parameters level
						# mean of MVN can be negative due to PCA - mean hprior
						hprior <- dunif(mmvn, min=-10000, max=10000, log=TRUE)
			
						if (model != "BM") { #OU hpriors
								  #here the entire vector is built (mean MVN, sig.sq MVN, alpha
								  hprior <- c(hprior, dgamma(c(smvn,bmou), shape=1.1, 
											 scale=5, log=TRUE)) 

						} else {# BM hpriors
								  #NULL is for prior on mu, which is not sampled
									
								 hprior <- c(hprior, dgamma(c(smvn, bmou), 
											 shape=1.1, scale=5, log=TRUE))  
											 
						}
						
						if (real.iter > 1) {
							 Lik <- LikA
							 priorMVN <- priorMVNA
							 priorBMOU <- priorBMOUA
						}
						
						# do this for first step always (because we need to have all probabiliiles)     
						if (r < update.freq[1] || real.iter == 1) { 
							 Lik <- likMultinorm(msp, ssp, traits, counts) # traits, counts
							 priorMVN <- likMVN(c(mmvn,smvn), msp, tree) # tree Conditional prior level
							 
							 if (model != "BM"){
								  priorBMOU <- likOU(bmou, ssp, tree, map) #  tree, map
							 } else {
								  priorBMOU <- likBM(bmou, ssp, tree) # tree
							 } 
							 
						} else if (r<update.freq[2]) {
							 priorMVN <- likMVN(c(mmvn, smvn), msp, tree)
						} else {
							 if (model != "BM") {
								  priorBMOU <- likOU(bmou, ssp, tree, map) #  tree, map
							 } else {
								  priorBMOU <- likBM(bmou, ssp, tree) # tree
							 } # Likelihood level
						}

						# Posterior calculation

						# jsut for 1 real.iter we need to copy all calculated likelihoods and priors
						if (real.iter == 1) {
							 LikA <- Lik
							 priorMVNA <- priorMVN
							 priorBMOUA <- priorBMOU
							 hpriorA <- hprior
							 postA <- (sum(Lik) + priorMVN + priorBMOU * temperature + sum(hprior))
						}
						post <- (sum(Lik) + priorMVN + priorBMOU * temperature + sum(hprior))
						# acceptance probability
						tryCatch(
						{
						if (post - postA + hasting.ratio >= log(runif(1))){
							 acc = acc + 1
							 LikA = Lik
							 priorMVNA = priorMVN
							 priorBMOUA = priorBMOU
							 hpriorA = hprior
							 postA = post
							 mspA = msp
							 sspA = ssp
							 mmvnA = mmvn
							 smvnA = smvn
							 bmouA = bmou
							}
						}
						,error = function(e) NULL
						)
						# log to file with frequency sampling.freq
						if (real.iter == 1){
							cat("generation",'\t',"posterior",'\n')
							if (model == "BM"){
							  cat(sprintf("%s\t", c("real.iter", "postA", "log.lik", 
										  "priorMBM", "priorVBM",  "sumHpriorA", "mbm_anc.st", 
										  "mbm_sig.sq", "vbm_sig.sq", "vbm_anc.st", 
										  paste("sp",seq(1:length(counts)),"_mean",sep=""), 
										  paste("sp",seq(1:length(counts)),"_var",sep=""),
										  "acc", "temperature")), "\n", append=FALSE, file=log.file)
							} else {
							  cat(sprintf("%s\t", c("real.iter", "postA", "log.lik", 
										  "priorMBM", "priorVOU",  "sumHpriorA", "mbm_anc.st", 
										  "mbm_sig.sq", "vou_alpha", "vou_sig.sq", "vou_anc.st", 
										  paste("vou_theta", seq(1:(length(bmouA)-3)),sep=""),
										  paste("sp",seq(1:length(counts)),"_mean",sep=""),
										  paste("sp",seq(1:length(counts)),"_var",sep=""),
										  "acc", "temperature")), "\n", append=FALSE, file=log.file)
							} 
							 
						}
						
						if (real.iter %% sampling.freq == 0 & real.iter >= burnin) {
							 cat(sprintf("%s\t", c(real.iter, postA, sum(LikA),
							 priorMVNA, priorBMOUA, sum(hpriorA), mmvnA, smvnA,
							 bmouA, mspA, sspA, (acc/iteration), temperature)),
							 "\n", append=TRUE, file=log.file) 
						}
						
						if (real.iter %% print.freq == 0 & real.iter >= burnin) {
							cat(real.iter,'\t',postA,'\n') 
						}
						  
				   real.iter = real.iter + 1
				   }
		 }# end of temperature
}



