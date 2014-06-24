
load("C:\\thematic\\HBM\\simulations\\!main_results\\params\\OU1\\rda\\sim107.rda")


load("C:\\thematic\\HBM\\simulations\\!main_results\\params\\OUM\\rda\\sim105.rda")

setwd("c:\\thematic\\github\\jive\\source\\")


load("C:\\thematic\\HBM\\simulations\\!main_results\\params\\BM\\rda\\sim103.rda")

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

# how often obtain
# temporary tesing for MCMC
jive.obj<-list()
jive.obj$traits<-res$data
jive.obj$counts<-res$observ.counts
jive.obj$tree<-res$phy
jive.obj$vcv<-vcv(res$phy)
jive.obj$map<-res$simmap
jive.obj$model<-res$model
jive.obj$nreg<-dim(res$simmap)[2]
str(jive.obj)


# ------------- BM
source("jive_mcmc_dev.R")
phy1 <- jive.obj$tree
tra1 <- jive.obj$traits
rownames(tra1) <- phy1$tip.label

my.jive <- make.jive(phy1, tra1,  model_var="BM", model_mean="BM")
jiveMCMC(my.jive, log.file="jive_mcmc_BM_xxx.log")

# ------------- OU1
source("jive_mcmc_dev.R")
phy1 <- jive.obj$tree
tra1 <- jive.obj$traits
rownames(tra1) <- phy1$tip.label

my.jive <- make.jive(phy1, tra1,  model_var="OU1", model_mean="BM")
jiveMCMC(my.jive, log.file="jive_mcmc_OU1.log")

# ------------- OUM
source("jive_mcmc_dev.R")
phy1 <- jive.obj$tree
tra1 <- jive.obj$traits
rownames(tra1) <- phy1$tip.label

n = length(phy1$tip.label)
states=c(rep(1,n%/%2),rep(2,(n - n%/%2)))
names(states)=phy1$tip.label
sim=make.simmap(phy1,states, model="SYM", nsim=1) # OUM model: output will be somehow random with same params due to states + simmap. 
my.jive <- make.jive(sim, tra1,  model_var="OUM", model_mean="BM")
jiveMCMC(my.jive, log.file="jive_mcmc_OUM.log")



ll[["traits"]]<-td$traits
ll[["counts"]]<-apply(td$traits, 1, function (x) {sum( !is.na(x) )})
ll[["tree"]]<-rel.sim(tr$tree)
ll[["model"]]<-model
ll[["nreg"]]<-dim(tr$tree)[2]



# make global variables from loaded rda.file





## run MCMC - potentially some of the MCMC params can be sent to bash script ?
#MCMC(update.freq=c(0.4,0.5),gbm_ou.ws=c(60,60,60,60,60), log_file=log.file, sampling_freq=1000, IT=2500000, burnin=0) 
n = 15
tree <- pbtree(b=1, n=n, scale=1, nsim=1, ape=TRUE)
map	 <- c(rep(1,n/2),rep(2,1+n/2))
names(map) <- tree$tip.label
simmap <- make.simmap(tree,map,model="SYM")
spec.obs=rpois(n, mean.obs)
M.spec.obs=max(spec.obs)
data.mat=matrix(rnorm(M.spec.obs*n,mean=mean.val,sd=sqrt(sigma.val)),nrow=n,ncol=M.spec.obs)
data.mat=cbind(as.matrix(M.spec.obs-spec.obs),data.mat)
foo<-function(x){
	to<-x[1]
	x[1:(to+1)]<-NA
	return(x[-1])
}

data.mat=t(apply(data.mat,1, foo))






res$data<-traits
res$observ.counts<-apply(res$data,1, function (x) {sum( !is.na(x) )})
	
	
SP_DATA = res$data #dataset - lines are vector of observations for each species
N_OBSERV = res$observ.counts #count of observations
TREE = res$phy
VCV.M = vcv(res$phy) #vcv matrix
MAP=res$simmap
MODEL=res$model 
NREG=dim(MAP)[2]

jive