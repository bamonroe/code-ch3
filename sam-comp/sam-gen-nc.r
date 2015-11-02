# Recreate Error Calculation
# Clear all
rm(list = ls())

#Libraries
library(microbenchmark)
library(Rcpp)
sourceCpp("~/code/Rcpp/halton.cpp")
sourceCpp("~/code/Rcpp/GRDD.cpp")
#sourceCpp("../Rcpp/DD.cpp")
#halton(init,H,prime)

library(parallel)
cores <- detectCores() / 2 
cores <- detectCores() -1
cores <- detectCores() 

GetRes <- function(sim){

	# The passed r value is the real "r", but all the calculations use 1-r
	# whenever r is needed, so subtract it once initially and save the 
	# computation times of all the subsequent calls, these tiny fractions
	# of a second add up with a very large snum.
	r1  <- 1 - sim[1]
	mu <- sim[2]

	# Calculate the 'Context'
	CTX <-  Max^r1/r1 - Min^r1/r1

	# Calculate the basline utility
	UA <- (pA0 * A0^r1/r1 + pA1 * A1^r1/r1)
	UB <- (pB0 * B0^r1/r1 + pB1 * B1^r1/r1)

	# Figure out which of these are errors
	A.err <- ifelse(UB > UA,1,0)
	B.err <- ifelse(UA > UB,1,0)

	# Make some vectors of certainty equivalents
	CEA <- (UA * r1)^(1/r1)
	CEB <- (UB * r1)^(1/r1)
	CEM <- ifelse(CEA>CEB,CEA,CEB)

	# Calculte the probability stuff, need to rebase the utility because 
	# of precision and machine limit problems. exp(>709) evalutes as Inf.
	
	UB.1 <- UB/CTX/mu - UA/CTX/mu
	UA.1 <- 0

	# Have things gone haywire because of the computer's inability to handle
	# numbers bigger than ~3e310 ?
	c.N <- is.nan(UB.1) | is.infinite(UB.1)

	pA <- ifelse( c.N ,		# Are we dealing with an insane number?
		# yes:
		ifelse( UB > UA , 0 , 1 ) ,
		# no, but are we making an insane number via exp?
		ifelse( UB.1 > 709 , 0 , exp(UA.1) / (1 + exp(UB.1)) )
	)

	pB <- 1 - pA

    # Gather the results 

	c(A.err , B.err , CEA , CEB , CEM , pA , pB)
    
}

DD <- function(cset, M,Errors,Cert,CEMax,Probs) {

	# c(A.err , B.err , CEA , CEB , CEM , pA , pB)

	EE <- sum(ifelse(cset == 0 , M[1:10], M[11:20]))

	PC <- prod(ifelse(cset == 0 , M[51:60], M[61:70]))

	WC <- mean(ifelse(cset == 0 , M[21:30] - M[31:40], M[31:40] - M[21:30]))

	WP <- mean(ifelse(cset == 0 , M[21:30] / M[41:50] , M[31:40] / M[41:50] ))

	# Errors[1:10,] is A
	# Errors[11:20,] is B

	# Make a vector that is 20 elements long, is equal to 1 if the choice is for that option
	# First 10 elements of Errors correspond to option A, second 10 to option B
	choice <- ifelse(cset==0,1,0)
	choice <- c(choice,cset)

	E.set <- apply(Errors,2,function(x) choice * x )
	E.num <- colSums(E.set)

	E.0 <- sum(E.num == 0) 
	E.1 <- sum(E.num == 1) 
	E.2 <- sum(E.num == 2) 
	E.3 <- sum(E.num == 3) 
	E.4 <- sum(E.num == 4) 
	E.5 <- sum(E.num == 5) 
	E.6 <- sum(E.num == 6) 
	E.7 <- sum(E.num == 7) 
	E.8 <- sum(E.num == 8) 
	E.9 <- sum(E.num == 9) 
	E.10 <- sum(E.num == 10) 

	#Cert   <- Res[21:40, ]
	#CEMax  <- rbind( Res[41:50, ], Res[41:50, ] )

#	WP.set <- apply(Cert,2,function(x) choice * x )
#	WP.set <- WP.set / CEMax
#	WP.num <- colMeans(WP.set) * 2
#
#	WP.1 <-	sum(ifelse(E.num == 1 , WP.num , 0 )) / E.1
#	WP.2 <-	sum(ifelse(E.num == 2 , WP.num , 0 )) / E.2
#	WP.3 <-	sum(ifelse(E.num == 3 , WP.num , 0 )) / E.3
#	WP.4 <-	sum(ifelse(E.num == 4 , WP.num , 0 )) / E.4
#	WP.5 <-	sum(ifelse(E.num == 5 , WP.num , 0 )) / E.5
#	WP.6 <-	sum(ifelse(E.num == 6 , WP.num , 0 )) / E.6
#	WP.7 <-	sum(ifelse(E.num == 7 , WP.num , 0 )) / E.7
#	WP.8 <-	sum(ifelse(E.num == 8 , WP.num , 0 )) / E.8
#	WP.9 <-	sum(ifelse(E.num == 9 , WP.num , 0 )) / E.9
#	WP.10 <- sum(ifelse(E.num == 10 , WP.num , 0 )) / E.10

	E.0 <- E.0 /snum 
	E.1 <- E.1 /snum
	E.2 <- E.2 /snum
	E.3 <- E.3 /snum
	E.4 <- E.4 /snum
	E.5 <- E.5 /snum
	E.6 <- E.6 /snum
	E.7 <- E.7 /snum
	E.8 <- E.8 /snum
	E.9 <- E.9 /snum
	E.10 <-E.10/snum 

	#data.frame(pattern=toString(cset),EE=EE,PC=PC,WC=WC,WP=WP,
	data.frame(EE = EE, PC = PC, WC = WC, WP = WP,
				E.0 = E.0,   
				E.1 = E.1,
				E.2 = E.2,
				E.3 = E.3,
				E.4 = E.4,
				E.5 = E.5,
				E.6 = E.6,
				E.7 = E.7,
				E.8 = E.8,
				E.9 = E.9,
				E.10 = E.10
	)

	#			WP.1 = WP.1,
	#			WP.2 = WP.2,
	#			WP.3 = WP.3,
	#			WP.4 = WP.4,
	#			WP.5 = WP.5,
	#			WP.6 = WP.6,
	#			WP.7 = WP.7,
	#			WP.8 = WP.8,
	#			WP.9 = WP.9,
	#			WP.10 = WP.10
}

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

	patlist <- t(pattern)

	D <- data.frame(DDcpp(pat=patlist, M=M,Errors=Errors,Cert=Cert,CEMax=CEMax,Probs=Probs ))
	#D <- lapply(patlist,DD, M=M,Errors=Errors,Cert=Cert,CEMax=CEMax,Probs=Probs)
	#D <- do.call(rbind,D)

	rm(Errors,Cert,CEMax,Probs)

	D[,5:15] <- D[,5:15] / snum

	colnames(D) <- c("EE","PC","WC","WP","E.0", "E.1", "E.2", "E.3", "E.4", "E.5", "E.6", "E.7", "E.8", "E.9", "E.10")

	D[,3:4] <- D[,3:4] / 10

#D$pattern <- apply(patlist,2,toString)

	D$rmean  <- rm
	D$mumean <- um
	D$rstd   <- rs
	D$mustd  <- us

	#D <- D[order(D$PC, decreasing=T),]

	M.EE  <- sum(D$PC*D$EE)
	M.WP  <- sum(D$PC*D$WP)
	M.WC  <- sum(D$PC*D$WC) 
	M.WE0 <- sum(D$PC*D$E.0)

	V.EE  <- sum(D$PC*(D$EE - M.EE )^2)
	V.WP  <- sum(D$PC*(D$WP - M.WP )^2)
	V.WC  <- sum(D$PC*(D$WC - M.WC )^2)
	V.WE0 <- sum(D$PC*(D$E.0 - M.WE0)^2)

	c(  M.EE=M.EE,
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

# The number of samples to draw
S <- 1000000

# Some boundary conditions
rm.min <- -5
rm.max <- 5
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

# Get rid of firs bunch of numbers to help combat possible correlations
burn <- 29
start <- 1 + burn
end <- S + burn

# Theta is our parameter matric
Theta <- data.frame(matrix(c(
	rm = halton(start,end,3)*rmdiff + rm.min,	# for r
	rs = halton(start,end,7)*rsdiff + rs.min,	
	um = halton(start,end,11)*umdiff + um.min,	# for mu 
	us = halton(start,end,13)*usdiff + us.min
	
	), 
	nrow=4,
	byrow=T
))

#Theta[2,] <- .15
#Theta[3,] <- .35
#Theta[4,] <- .15

# What is the number of simulations to run per sample
snum <- 10000

# Set up the Halton sequences
burn <- 30
start <- 1 + burn
end <- snum + burn

H1 <- halton(start,end,3)	# for r
H2 <- halton(start,end,7)	# for mu

# How many choices per subject
cnum <- length(A0)

# Define some population parameters, mean and standard devitions of r
# and fechner values

#microbenchmark(
#	MM  <- mclapply( Theta, getAgg,mc.cores=cores), 
#	times=20
#	)


MM  <- mclapply( Theta, getAgg,mc.cores=cores)
#MM  <- lapply( Theta, getAgg)
MM <- do.call(rbind,MM)
MM <- data.frame(MM)
row.names(MM) <- 1:nrow(MM)

plot(MM$rm,MM$M.WP)


aggdir <- "../data/aggdat/"
fname <- paste(aggdir,"Agg",S,"-sim",snum,"-",rm.min,rm.max,".Rda",sep="")
#fname <- paste("Test.Rda")
save(MM,file=fname)
print(fname)

#source("/home/woodape/Dropbox/Ongoing Projects/Thesis/Code/sam-comp/comp.r", echo = TRUE, print.eval = TRUE, chdir = TRUE)

