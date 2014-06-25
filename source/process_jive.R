require(TeachingDemos)
require(hdrcde)
require(akima)
require(coda)

setwd("c:\\thematic\\github\\jive\\source\\resulting_files\\")

jive.mode<-function(x){
	y=hdr(x, all.modes=F)$mode
	return(y)#print(y)
}

proc.jive<-function(log.file = "jive_mcmc_OU1.log", n = 50, stat=jive.mode, burning = 0, probHPD = 0.95, ...){
		
		options(scipen=999)
		ll 	= list()
		l	= read.csv(log.file, stringsAsFactors=FALSE, header=T, sep="\t")
		g	= dim(l)[1] # number of generations
		p   = dim(l)[2] # number of parameters
		
		if (burning == 0) {
			b = c(p%/%2, p) # keep only half of chain if no burning is given
		} else {
			b = c(burning, g)
		}
		
		# e.g. start of params(7)# end of params (12)# start of spsp params (13)# end of spsp params(112)		
		s = c(7, (p - 3 - n * 2), (p - 2 - n * 2), (p - 3)) 

		
		# apply function to calculate summary
		ll$prior_pars  = apply(l[b[1]:b[2], c(s[1]:s[2])], 2, stat) # mode for prior level pars
		ll$mean_pars    = apply(l[b[1]:b[2], c(s[3]:(s[3] + n - 1))], 2, mean) # mean for lik level pars (mode is too slow)
		ll$var_pars    = apply(l[b[1]:b[2], c((s[3] + n):s[4])], 2, mean)
		names(ll$mean_pars) = sub("_mean","",names(ll$mean_pars))
		names(ll$var_pars) = sub("_var","",names(ll$var_pars))

		# calculate mcmc and hpd (coda package)
		my.mcmc        = mcmc(data=l[b[1]:b[2], 1:(p-1)])
		my.hpd         = data.frame(HPDinterval(my.mcmc), prob = probHPD, ...)
		
		# get the HPD intervals
		ll$prior_hpd.l = data.frame(my.hpd)[c(s[1]:s[2]), 1]
		ll$prior_hpd.u = data.frame(my.hpd)[c(s[1]:s[2]), 2]
		
		ll$mean_hpd.l   = data.frame(my.hpd)[c(s[3]:(s[3] + n - 1)), 1]
		ll$mean_hpd.u   = data.frame(my.hpd)[c(s[3]:(s[3] + n - 1)), 2]
		
		ll$var_hpd.l   = data.frame(my.hpd)[c((s[3] + n):s[4]), 1]
		ll$var_hpd.u   = data.frame(my.hpd)[c((s[3] + n):s[4]), 2]
		
		names(ll$prior_hpd.l) <- names(ll$prior_pars)
		names(ll$prior_hpd.u) <- names(ll$prior_pars)	

		names(ll$mean_hpd.l) <- names(ll$mean_pars)
		names(ll$mean_hpd.u) <- names(ll$mean_pars)			
		
		names(ll$var_hpd.l) <- names(ll$var_pars)
		names(ll$var_hpd.u) <- names(ll$var_pars)	
		
		return(ll)


}

plot.jive <- function(tree, proc.jive, regime, cols=c("blue", "green"), cex.label=0.7, cex.circle=2, lab.off = 0.025, ladder=TRUE, ...){

		if (ladder){
			tree = ladderize(tree, right = TRUE)
		}
		regime <- regime[match(tree$tip.label, names(regime))]
		traits  <- proc.jive$var_pars[match(tree$tip.label, names(proc.jive$var_pars))]
		
		plot(tree,  label.offset = lab.off, show.tip.label = TRUE, cex = cex.label)
		tiplabels(pch = 19, col = cols[regime], cex = cex.circle * (traits/max(traits)))

}

plot.jive(phy1, my.l, regime)

foo <- function(x){
	#plot(density(x[1],x[2]), axes = FALSE, xlab=NULL, ylab=NULL, main=NULL)
	plot(density(x[1],x[2]), main="", xlab="", ylab="", xlim=c(340,360), xaxt="n", yaxt="n", axes=F)
	abline(v=350)

}

par(mfrow=c(5,5))
apply(dd, 1, foo)

