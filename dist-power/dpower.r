# Doing some power calculations for sampling from the distributions generated
# in Agg-Gen

# Clear all
rm(list=ls())

# Libraries
library(MASS)
library(parallel)

cores <- detectCores() - 1

# Load in Agg-Data
load("../data/agg-dat/Agg30k-S500.Rda")

# Get rid of everything but distributions
A <- data.frame(matrix(c(rm=MM$rm,rs=MM$rs,um=MM$um,us=MM$us),byrow=TRUE,nrow=4))

rm(MM)


#k.m <- mean( (A[3,]^2)/(A[4,]^2) )
#t.m <- mean( (A[4,]^2)/A[3,] )

getNSam <- function(theta,n){

	rm <- theta[1]
	rs <- theta[2]
	um <- theta[3]
	us <- theta[4]

	# Fechner will use a gamma distribution, so need to back out shape and
	# scale parameters
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	f <- matrix(c(
			r = rnorm(n=n,mean=rm,sd=rs),
			u = rgamma(n=n,shape=k,scale=t)
		),ncol=2)


	r <- rnorm(n=n,mean=rm,sd=rs)

	fit <- fitdistr(f[,1],"normal")

	res <- c(m.est=fit$estimate[1], s.est=fit$estimate[2], m.se=fit$sd[1], s.se=fit$sd[2], rm=rm, rs=rs)

	return(res)
    
}

getAgg <- function(n){
    
	fits  <- lapply(A,getNSam,n=n)
#	fits  <- mclapply(A,getNSam,n=n,mc.cores=cores)

	full  <- data.frame(do.call(rbind,fits))

	full$rm.diff <- abs(full$rm - full$m.est.mean)
	full$rs.diff <- abs(full$rs - full$s.est.sd)

	den.m <- mean(full$rm.diff)
	den.s <- mean(full$rs.diff)


	dis.m <- fitdistr(full$rm.diff,"gamma")
	dis.s <- fitdistr(full$rs.diff,"gamma")

	res <- c(n,dis.m$estimate[1], dis.s$estimate[2], dis.m$sd[1], dis.s$sd[2])

	return(res)

}

sam.size <- seq(from=10,to=500,by=10)

agg.sams <- mclapply(sam.size,getAgg,mc.cores=cores)
#agg.sams <- lapply(sam.size,getAgg)

err <- data.frame(do.call(rbind,agg.sams))
colnames(err) <- c("N","shape", "rate","shape.se","rate.se")

err$scale <-  1 / err$rate

err$mean <- err$shape * err$scale
err$sd   <- (err$shape * err$scale^2 ) ^.5

