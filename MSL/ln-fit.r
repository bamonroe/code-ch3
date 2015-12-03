# There is no analytical solution to back out the parameters for the normal
# distribution from the logit-normal distribution. This is an attempt to get
# them from numberical optimization.

#libraries
library(Rcpp)

# Source the cpp stuff
sourceCpp("../Rcpp/halton.cpp")

# Get the data
load("../data/choice-dat/choice10.Rda")

rm <- mean(D$r)
rs <- sd(D$r)

um <- mean(D$mu)
us <- sd(D$mu)

k <- um^2 / us^2
t <- us^2 / um

lu <- function(par){
	lum <- par[1]
	lus <- exp(par[2])
	lust <- exp(par[3])

	oo <- exp(qnorm(halton(1,500,3),lum,lus,lust))
	oo <- oo / (1 + oo)
	oo * lust
}

lr <- function(par){
	lrm <- par[1]
	lrs <- exp(par[2])
	lrst <- exp(par[3])
	lrsh <- par[4]

	oo <- exp(qnorm(halton(1,500,3),lrm,lrs))
	oo <- oo / (1 + oo)
	oo * lrst + lrsh
}

tr <- function(par){
	lrm <- par[1]
	lrs <- exp(par[2])
	lrst <- exp(par[3])
	lrsh <- par[4]

	c(lrm,lrs,lrst,lrsh)
} 

tu <- function(par){
	lum <- par[1]
	lus <- exp(par[2])
	lust <- exp(par[3])

	c(lum,lus,lust)

}

mingamma <- function(par){

	lum <- par[1]
	lus <- exp(par[2])
	lust <- exp(par[3])

	gamma <- qgamma(halton(30,531,3),k,1/t)

	lgamma <- exp(qnorm(halton(30,531,3),lum,lus))
	lgamma <- lgamma / (1 + lgamma)
	lgamma <- lgamma * lust

	gamma  <- gamma[order(gamma)]
	lgamma <- lgamma[order(lgamma)]

	sum((gamma - lgamma)^2)
}

minnorm <- function(par){

	lrm <- par[1]
	lrs <- exp(par[2])
	lrst <- exp(par[3])
	lrsh <- par[4]

	norm <- qnorm(halton(30,531,3),rm,rs)

	lnorm <- exp(qnorm(halton(30,531,3),lrm,lrs))
	lnorm <- lnorm / (1 + lnorm)
	lnorm <- lnorm * lrst + lrsh

	norm  <- norm[order(norm)]
	lnorm <- lnorm[order(lnorm)]

	sum((norm - lnorm)^2)
}

dist <- c(-2,.5,15,0,-4,.6,14)
init <- c(dist[1],log(dist[2]),log(dist[3]),dist[4],dist[5],log(dist[6]),log(dist[7]))

init.n <- init[1:4]

res <- optim(init.n,minnorm)

pr <- tr(res$par)
lnorm <- lr(res$par)
plot(density(lnorm))

print(c(mean(lnorm),sd(lnorm)))


init.g <- init[5:7]

res <- optim(init.g,mingamma)

pq <- tu(res$par)
lgamma <- lu(res$par)
plot(density(lgamma))

print(c(mean(lgamma),sd(lgamma)))

