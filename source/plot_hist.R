# get true dist
d  <- rbind(rnorm(1000, mean=6, sd=0.01), rnorm(1000, mean=3, sd=1.5))

# rescale them
rs <- rescale(d, to=c(10.5,11.5))
rs1 <- re[1:1000]
rs2 <- re[10001:2000]

# get density coord
de1 <- density(rs1)
xx1 <- de1$xx
yy1 <- de1$yy

de2 <- density(rs2)
xx2 <- de2$xx
yy2 <- de2$yy


layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))


plot(t1, show.tip.label = FALSE)
getTreeCoords <- function(tree, xwidth=1, yheight=0.4){
	
	d = list()
	h = vcv(tree)[1] # height
	n = dim(vcv(tree))[2] # num spec

	d$xxl <- rep(h, n) + 0.5
	d$xxu <- rep(h+xwidth, n)
	
	d$yyl <- 1:n # - yheight
	d$yyu <- d$yyl + 2*yheight
	
	return(d)
}

getCoordNDC <- function(d){

	d$xxl<-grconvertX(d$xxl, from = "user", to = "ndc")
	d$xxu<-grconvertX(d$xxu, from = "user", to = "ndc")
	
	d$yyl<-grconvertX(d$yyl, from = "user", to = "ndc")
	d$yyu<-grconvertX(d$yyu, from = "user", to = "ndc")
	
	return(d)
}
rr1 <- getTreeCoords(t1)
rr <- getCoordNDC(rr1)


#par(fig=c(rr$xxl[i],rr$xxu[i], rr$yyl[i], rr$yyu[i]),new = TRUE, mar = rep(0, 4))
#plot(0:1, 0:1, las = 1, main="", xlab="", ylab="", xaxt="n", yaxt="n", axes=F )

for (i in 1:10){
	par(fig=c(rr$xxl[i],rr$xxu[i], rr$yyl[i], rr$yyu[i]), new=TRUE, mar = rep(0, 4))
	#plot(0:1, 0:1, las = 1, main="", xlab="", ylab="", xaxt="n", yaxt="n", axes=F )
	hist(rnorm(1000,0.5,0.1), main="", xlab="", ylab="", xaxt="n", yaxt="n", axes=F, xlim=c(0,1))
	
}


require(TeachingDemos)

plot(t1, show.tip.label = FALSE)
getTreeCoords <- function(tree, xwidth=1, yheight=0.4){
	
	d = list()
	h = vcv(tree)[1] # height
	n = dim(vcv(tree))[2] # num spec

	d$xxl <- rep(h, n) + 0.5
	d$xxu <- rep(h+xwidth, n)
	
	d$yyl <- 1:n # - yheight
	d$yyu <- d$yyl + 2*yheight
	
	return(d)
}
rr1 <- getTreeCoords(t1)

for (i in 1:10){
	if (i < 4){
		par(mgp=c(0, 0, 0))
		subplot( hist(rnorm(100, sd=0.1),xlab='',ylab='',main='', yaxt="n", xlim=c(-2,2), tck=-0.05, cex.axis=0.5), rr1$xxl[i], rr1$yyl[i], size=c(0.5,0.5))
		
	}
	else{
		subplot( hist(rnorm(100, sd=1),xlab='',ylab='',main='', yaxt="n",  xlim=c(-2,2), tck=-0.05, cex.axis=0.5), rr1$xxl[i], rr1$yyl[i], size=c(0.5,0.5))
	}
}


subplot( hist(rnorm(100),xlab='',ylab='',main='', xaxt="n", yaxt="n", axes=F), rr1$xxl[i], rr1$yyl[i], size=c(0.5,0.5))


