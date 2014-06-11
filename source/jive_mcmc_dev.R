### MCMC TEMPLATE ###

# LOAD PACKAGES
library(OUwie)
library(nloptr)
library(ape)
library(MASS)
library(phytools)


source("lik_bm.R")
source("lik_ou.R")
source("lik_multinorm.R")
source("hpriors.R")
source("proposals.R")
source("utils.R")
source("jive_prep_dev.R")


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

					




# MCMC MAIN ALGORITHM
jiveMCMC <- function(jive, log.file="jive_mcmc.log", sampling.freq=1000, print.freq=100, 
				ncat=1, beta.param=0.3, ngen=5500000, burnin=0, update.freq=NULL, model="OU")
{
	
			 
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
							
							# initialize MCMC parameters
							mspA  <- jive$lik$mspinit # initialize means for species
							sspA  <- jive$lik$sspinit # initialize sigma.sq for species
							mmvnA <- jive$prior_mean$init[1] # could be a random number, initialize mean for MVN
							smvnA <- jive$prior_mean$init[2] # could be a random number, initialize sigma.sq for MVN
							bmouA <- jive$prior_var$init # could be aither more realistic values such as means and sds of true data
							 
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
								  msp[ind] <- slidingWin(jive$lik$mspinit[ind], jive$lik$mspws[ind]) 
							 } else {
								  # updating random 5 values from the vector of sigmas
								  ssp[ind] <- abs(slidingWin(jive$lik$sspinit[ind], jive$lik$sspws[ind]))
							 }
						}
						# update MVN parameters
						else if (r < update.freq[2]) {
						
							 if (runif(1) < 0.5){
								  # updating mmvn - negative values allowed
								  mmvn <- slidingWin(jive$prior_mean$init[1], jive$prior_mean$ws[1]) 
							 } else {
								  # updating smvn
								  smvn <- abs(slidingWin(jive$prior_mean$init[2], jive$prior_mean$ws[2]) )
							 }
							 
						} else {# update BMOU parameters 
												 
							 ind <- sample(1:length(jive$prior_var$ws), 1)
							 bmou[ind] <- abs(slidingWin(jive$prior_var$init[ind], jive$prior_var$ws[ind])) # updating bmou parameters
							 
						}
						

						# Hyperpriors on all parameters level
						# mean of MVN can be negative due to PCA - mean hprior
											
						hprior_mean <- mapply(do.call, jive$prior_mean$hprior, lapply(c(mmvn, smvn), list))
						hprior_var <- mapply(do.call, jive$prior_var$hprior, lapply(bmou, list))
						hprior <- c(hprior_mean, hprior_var)
						#print(hprior)
						
						if (real.iter > 1) {
							 Lik <- LikA
							 Prior_mean <- Prior_meanA #priorMVNA
							 Prior_var <- Prior_varA #priorBMOUA
						}
						
						# do this for first step always (because we need to have all probabiliiles)     
						if (r < update.freq[1] || real.iter == 1) {
						
							 Lik 		<- jive$lik$model(msp, ssp, jive$data$traits, jive$data$counts) # traits, counts
							 #print(Lik)
							 Prior_mean <- jive$prior_mean$model(c(mmvn,smvn), msp, jive$data$tree) # tree Conditional prior level
							 #print(paste("Prior mean ", Prior_mean))
							 Prior_var 	<- jive$prior_var$model(bmou, ssp, jive$data$tree, jive$data$map)
							 #print(paste("Prior var ", Prior_var))
							 										 
						} else if (r<update.freq[2]) {
							 Prior_mean <- jive$prior_mean$model(c(mmvn, smvn), msp, jive$data$tree)
						} else {
							 Prior_var 	<- jive$prior_var$model(bmou, ssp, jive$data$tree, jive$data$map)
						}

						# Posterior calculation

						# jsut for 1 real.iter we need to copy all calculated likelihoods and priors
						if (real.iter == 1) {
							 LikA <- Lik
							 Prior_meanA <- Prior_mean
							 Prior_varA <- Prior_var
							 hpriorA <- hprior
							 postA <- (sum(Lik) + Prior_mean + Prior_var * temperature + sum(hprior))
						}
						post <- (sum(Lik) + Prior_mean + Prior_var * temperature + sum(hprior))
						# acceptance probability
						tryCatch(
						{
						if (post - postA + hasting.ratio >= log(runif(1))){
							 acc = acc + 1
							 LikA = Lik
							 Prior_meanA = Prior_mean
							 Prior_varA = Prior_var
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
										  "Prior_mean", "Prior_var",  "sumHpriorA", "mbm_anc.st", 
										  "mbm_sig.sq", "vbm_sig.sq", "vbm_anc.st", 
										  paste("sp",seq(1:length(jive$data$counts)),"_mean",sep=""), 
										  paste("sp",seq(1:length(jive$data$counts)),"_var",sep=""),
										  "acc", "temperature")), "\n", append=FALSE, file=log.file)
							} else {
							  cat(sprintf("%s\t", c("real.iter", "postA", "log.lik", 
										  "Prior_mean", "Prior_var",  "sumHpriorA", "mbm_anc.st", 
										  "mbm_sig.sq", "vou_alpha", "vou_sig.sq", "vou_anc.st", 
										  paste("vou_theta", seq(1:(length(bmouA)-3)),sep=""),
										  paste("sp",seq(1:length(jive$data$counts)),"_mean",sep=""),
										  paste("sp",seq(1:length(jive$data$counts)),"_var",sep=""),
										  "acc", "temperature")), "\n", append=FALSE, file=log.file)
							} 
							 
						}
						
						if (real.iter %% sampling.freq == 0 & real.iter >= burnin) {
							 cat(sprintf("%s\t", c(real.iter, postA, sum(LikA),
							 Prior_meanA, Prior_varA, sum(hpriorA), mmvnA, smvnA,
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



