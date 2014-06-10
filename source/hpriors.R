

make.hpfunction <-function(hpf="Uniform", hp.pars){


		if (hpf == "Uniform"){
			my.f <- function(x){
				hp <- sum(dunif(x, min=hp.pars[1], max=hp.pars[2], log=TRUE))
				return(hp)
			}
			
		}
		
		if (hpf == "Gamma"){
			my.f <- function(x){
				hp <- sum(dgamma(x, shape=hp.pars[1], scale=hp.pars[2], log=TRUE))
				return(hp)
			}
		}
		
		if (hpf == "Normal"){
			my.f <- function(x){
				hp <- sum(dnorm(x, mean=hp.pars[1], sd=hp.pars[2], log=TRUE))
				return(hp)
			}
		}	
		
		
	return(my.f)
}