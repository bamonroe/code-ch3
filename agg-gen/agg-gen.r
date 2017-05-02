# Recreate Error Calculation
# Clear all
rm(list = ls())

#Libraries
library(ctools)
c.library("halton","Rcpp")

# Compile C++ code
c.sourceCpp("../Rcpp/sam-gen.cpp")

small.run <- T

getAgg <- function(theta){

	rm <- theta[1]
	rs <- theta[2]
	um <- theta[3]
	us <- theta[4]

	# Fechner will use a gamma distribution, so need to back out shape and
	# scale parameters
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	sim <- matrix(  c(qnorm(H1,mean = rm, sd = rs), qgamma(H2,shape = k, scale = t) ) ,byrow=T, ncol=snum  )

	# Give the very wide bounds we've mad for ourselves, occasionally we get 
	# a fechner value equal to 0. This creates an "oh shi-" problem later on.
	# Keep replacing 0s with another random draw from the same distribution until
	# they go away.
	while(any(sim[,2]==0)){
		sim[,2] <- ifelse(sim[,2]==0 , rgamma(1,shape = k, scale = t),sim[,2] )
	}
	while(anyNA(sim[,2])){
		sim[,2] <- ifelse(is.na(sim[,2]) , rgamma(1,shape = k, scale = t),sim[,2] )
	}

	Res <- getRes(sim,A0,A1,B0,B1,
				  pA0,pA1,pB0,pB1,
				  Max,Min)

	#c(A.err , B.err , CEA , CEB , CEM , pA , pB)

	M <- rowMeans(Res)

	Errors <- Res[1:20, ]
	Cert   <- Res[21:40, ]
	CEMax  <- rbind(Res[41:50, ] , Res[41:50, ])
	Probs  <- Res[51:70, ]

	rm(Res,sim)

	# Returns a dataset with 16 columns.
	D <- data.frame(DDcpp(pat=patlist, M=M,Errors=Errors,Cert=Cert,CEMax=CEMax,Probs=Probs ))

	# The first five are
	#EE  -  Expected number of errors
	#PC  -  Simulated probability
	#LPC -  Log-Simulated probability !! Not the same as log(PC) !!
	#WC  -  Welfare Surplus metric
	#WP  -  Welfare proportion metirc

	# The next 10 are the total number of agents with [0-10] errors, so we need to divide it
	# by the number of agents we simulated to get the probability of [0-10] errors
	D[,6:16] <- D[,6:16] / snum

	colnames(D) <- c("EE","PC","LPC","WC","WP","E.0", "E.1", "E.2", "E.3", "E.4", "E.5", "E.6", "E.7", "E.8", "E.9", "E.10")

	#
	D[,4:5] <- D[,4:5] / 10

	M.LPC  <- sum(D$PC*D$LPC)
	M.EE  <- sum(D$PC*D$EE)
	M.WP  <- sum(D$PC*D$WP)
	M.WC  <- sum(D$PC*D$WC) 
	M.WE0 <- sum(D$PC*D$E.0)

	V.LPC  <- sum(D$PC*(D$LPC - M.LPC )^2)
	V.EE  <- sum(D$PC*(D$EE - M.EE )^2)
	V.WP  <- sum(D$PC*(D$WP - M.WP )^2)
	V.WC  <- sum(D$PC*(D$WC - M.WC )^2)
	V.WE0 <- sum(D$PC*(D$E.0 - M.WE0)^2)

	c(M.LPC = M.LPC,  
		V.LPC=V.LPC,
		M.EE=M.EE,
		V.EE=V.EE,
		M.WP=M.WP,
		V.WP=V.WP,
		M.WC=M.WC,
		V.WC=V.WC,
		M.WE0=M.WE0,
		V.WE0=V.WE0,
		rm=theta[1],
		rs=theta[2],
		um=theta[3],
		us=theta[4])

}

pattern <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
patlist <- t(pattern)

c.export("patlist")

# Set up our instrument
A0 <- rep(1.60,10)
A1 <- rep(2.00,10)
B0 <- rep(0.10,10)
B1 <- rep(3.85,10)

pA0 <- seq(.9,0,by=-.1)
pB0 <- seq(.9,0,by=-.1)
pA1 <- seq(.1,1,by=.1)
pB1 <- seq(.1,1,by=.1)

Max <- rep(3.85,10)
Min <- c(rep(0.10,9),2)

c.export("A0","A1","B0","B1","pA0","pA1","pB0","pB1","Max","Min")

# The number of samples to draw
S <- ifelse(small.run, 1000, 500000)

# Some boundary conditions
rm.min <- -1.9		# Just below the HL-MPL
rm.max <- 1.55		# Just above the HL-MPL
rs.min <- .01
rs.max <- 1
um.min <- .05
um.max <- 2.25
us.min <- .01
us.max <- .75

rm <- runif(S,min=rm.min,max=rm.max)
rs <- runif(S,min=rs.min,max=rs.max)
um <- runif(S,min=um.min,max=um.max)
us <- runif(S,min=us.min,max=us.max)

# Need to get the size of the interval to transform the (0,1) interval
# into a wider interval
rmdiff <- rm.max - rm.min
rsdiff <- rs.max - rs.min

umdiff <- um.max - um.min
usdiff <- us.max - us.min

# Theta is our parameter matric
Theta <- data.frame(matrix(c(
	rm = halton(S,3)*rmdiff + rm.min,	# for r
	rs = halton(S,7)*rsdiff + rs.min,	
	um = halton(S,11)*umdiff + um.min,	# for mu 
	us = halton(S,13)*usdiff + us.min
	
	), 
	nrow=4,
	byrow=T
))

# What is the number of simulations to run per sample
snum <- ifelse(small.run, 100, 10000)

# Set up the Halton sequences
H1 <- halton(snum,3)	# for r
H2 <- halton(snum,7)	# for mu

c.export("H1","H2")

# Define some population parameters, mean and standard devitions of r
# and fechner values

MM  <- c.lapply( Theta, getAgg)
#MM  <- lapply( Theta, getAgg)
MM <- do.call(rbind,MM)
MM <- data.frame(MM)
row.names(MM) <- 1:nrow(MM)

aggdir <- "../data/agg-dat/"
fname <- paste(aggdir,"Agg-sim.Rda",sep="")
#fname <- paste(aggdir,"Agg",S,"-sim",snum,"-",rm.min,rm.max,".Rda",sep="")
#fname <- paste("Test.Rda")
save(MM,file=fname)
print(fname)

