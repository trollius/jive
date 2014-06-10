#ll[["traits"]] <- td$data
#ll[["counts"]] <- apply(td$data, 1, function (x) {sum( !is.na(x) )})
#ll[["tree"]]   <- td$phy
#ll[["vcv"]]    <- vcv(td$phy)
#ll[["model"]]  <- model
#ll[["nreg"]]   <- dim(td$phy)[2]

class(ll) <- 'jive'


traits.jive <- function(jive.object)
{
	return(jive.object[['traits']])
}
 
traits <- function(jive)
{
	UseMethod('traits', jive)
}
 
model.jive  <- function(jive.object)
{
	return(jive.object[['model']])
}
 
model <- function(jive)
{
	UseMethod('model', jive)
}

model(my.jive)

"model<-"<-function(x,value){
	UseMethod("model<-")
}

"model<-.jive"<-function(x,value){
	x$model<-value
	x
}


model(my.jive) <- 'TTR'