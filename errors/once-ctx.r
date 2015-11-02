# Recreate Error Calculation
# Clear all
#rm(list = ls())

library(Rcpp)
sourceCpp("~/git/thesis/halton/halton.cpp")
#halton(init,H,prime)

library(parallel)
cores <- detectCores() - 1
cores <- 4

# Make the pattern
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

# What is the number of simulations to run
snum <- 1000000
#snum <- 100
# How many choices per subject
cnum <- length(A0)

start <- proc.time()

# For every simulation, I want the number of errors. I also want this per
# pattern, so every subject will get a 

# Define some population parameters, mean and standard devitions of r
# and fechner values

rm <- .65
rs <- .3
um <- .35
us <- .3

# Fechner will use a gamma distribution, so need to back out shape and
# scale parameters
k <- (um^2)/(us^2)
t <- (us^2)/um

burn <- 60
start <- 1 + burn
end <- snum + burn

H1 <- halton(start,end,3)
H2 <- halton(start,end,7)

sim <- data.frame( matrix(  c(qnorm(H1,mean = rm, sd = rs), qgamma(H2,shape = k, scale = t) ) ,byrow=T, ncol=snum  ) )

rm(H1,H2)

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

patlist <- data.frame(t(pattern))

GetRes <- function(sim,cset){

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
	c.AB <- UA > UB & cset == 1
	c.BA <- UB > UA & cset == 0
	e1 <- c.AB | c.BA

	c.err <- ifelse(e1,1,0)

	# Make some vectors of certainty equivalents
	CEA <- (UA * r1)^(1/r1)
	CEB <- (UB * r1)^(1/r1)

	# Get the CE of the chosen option
	W0 <- ifelse(cset==0,CEA,CEB)

	# Get th CE of the unchosen option
	W1 <- ifelse(cset==1,CEA,CEB)

	# Get the CE of the Error only
	WE <- ifelse(c.err==1,W0-W1,0)

	# Get the maximum CE
	WMax <- ifelse(CEA>CEB,CEA,CEB)

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
		ifelse( UB.1 > 709 , 0 , exp(UA.1) / (exp(UA.1) + exp(UB.1)) )
	)

	pB <- 1 - pA

	# Store the probabilities of the chosen options
	p0 <- ifelse(cset==0,pA,pB)

    # Gather the results 
    res <- c(0,0,0,0,0)

	# How many errors
    res[1] <- sum(c.err)

    # Probability of this pattern
    res[2] <- prod(p0)

    # Welfare metric mirroring opportunity cost of whole choice set
    res[3] <- sum(W0 - W1)

    # Welfare metric for proportion of welfare obtained / total availavle
	res[4] <- sum(W0)/sum(WMax)

	# Welfare metric opportunity cost of error only
	res[5] <- sum(WE)

	return(res)

	# 1- Err# , 2-Prob, 3-WOC, 4-Wprop
    
}



DD <- function(pat) {
    
	print(toString(pat))

	cset <- matrix(pat,nrow=10)

	# 1- Err# , 2-Prob, 3-WOC, 4-Wprop

	result <- lapply(sim,GetRes,cset=cset)
	result <- matrix( unlist(result), nrow=5 )

	EE <- sum(result[1,])/snum
	PC <- sum(result[2,])/snum

	WC <- sum(result[3,])/snum
	WP <- sum(result[4,])/snum
	WE <- sum(result[5,])/snum

	E.0 <- sum(result[1,] == 0)
	E.1 <- sum(result[1,] == 1)
	E.2 <- sum(result[1,] == 2)
	E.3 <- sum(result[1,] == 3)
	E.4 <- sum(result[1,] == 4)
	E.5 <- sum(result[1,] == 5) 
	E.6 <- sum(result[1,] == 6)
	E.7 <- sum(result[1,] == 7)
	E.8 <- sum(result[1,] == 8)
	E.9 <- sum(result[1,] == 9)
	E.10 <-sum(result[1,] == 10)

	WP.1 <-	sum(ifelse(result[1,] == 1, result[4,],0 )) / E.1
	WP.2 <-	sum(ifelse(result[1,] == 2, result[4,],0 )) / E.2
	WP.3 <-	sum(ifelse(result[1,] == 3, result[4,],0 )) / E.3
	WP.4 <-	sum(ifelse(result[1,] == 4, result[4,],0 )) / E.4
	WP.5 <-	sum(ifelse(result[1,] == 5, result[4,],0 )) / E.5
	WP.6 <-	sum(ifelse(result[1,] == 6, result[4,],0 )) / E.6
	WP.7 <-	sum(ifelse(result[1,] == 7, result[4,],0 )) / E.7
	WP.8 <-	sum(ifelse(result[1,] == 8, result[4,],0 )) / E.8
	WP.9 <-	sum(ifelse(result[1,] == 9, result[4,],0 )) / E.9
	WP.10 <- sum(ifelse(result[1,] == 10, result[4,],0 )) / E.10



	data.frame(pattern=toString(pat),
	  PC=PC,
	  WC=WC,
	  WP=WP,
	  WE=WE,
	  EE=EE,
	  E.0=E.0,
	  E.1=E.1,
	  E.2=E.2,
	  E.3=E.3,
	  E.4=E.4,
	  E.5=E.5,
	  E.6=E.6,
	  E.7=E.7,
	  E.8=E.8,
	  E.9=E.9,
	  E.10=E.10,
	  WP.1=WP.1,
	  WP.2=WP.2,
	  WP.3=WP.3,
	  WP.4=WP.4,
	  WP.5=WP.5,
	  WP.6=WP.6,
	  WP.7=WP.7,
	  WP.8=WP.8,
	  WP.9=WP.9,
	  WP.10=WP.10
	)

}

start <- proc.time()
	#D <- lapply(patlist,DD)
	D <- mclapply(patlist,DD,mc.cores=cores)
	D <- do.call(rbind,D)

end <- proc.time()

D$rmean  <- rm
D$mumean <- um
D$rstd   <- rs
D$mustd  <- us
D$ID     <- pattern[,1]

save(D ,file=paste("Dat-CTX.Rda",sep=""))

