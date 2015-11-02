# Recreate Error Calculation
# Clear all
rm(list = ls())

#Libraries
library(microbenchmark)
library(Rcpp)
sourceCpp("../Rcpp/halton.cpp")
#halton(init,H,prime)
sourceCpp("../Rcpp/GRDD.cpp")

library(parallel)
cores <- detectCores() / 2 
cores <- detectCores() 

# Pattern matrix, where each column is a choice pattern
patmat <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
patmat <- t(patmat)

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

# How many simulations to run per sample
snum <- 2000000

# Set up the Halton sequences
burn <- 30
start <- 1 + burn
end <- snum + burn

H1 <- halton(start,end,3)	# for r
H2 <- halton(start,end,7)	# for mu

# How many choices per subject
cnum <- length(A0)

# Define the distributional parameters
rm <- .65
rs <- .3
um <- .35
us <- .3

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

# getRes is written in C++ and sourced from GRDD.cpp
Res <- getRes(sim,A0,A1,B0,B1,
				pA0,pA1,pB0,pB1,
				Max,Min)

# Free up some memory
rm(sim)

# Each column of Res is created from the following combined vector
#c(A.err , B.err , CEA , CEB , CEM , pA , pB)

M <- rowMeans(Res)

Errors <- Res[1:20, ]
Cert   <- Res[21:40, ]
CEMax  <- rbind(Res[41:50, ] , Res[41:50, ])
Probs  <- Res[51:70, ]

# Free up some memory
rm(Res)

# DDcpp is written in C++ and sourced from GRDD.cpp
D <- data.frame( DDcpp(pat=patmat,M=M,Errors=Errors,Cert=Cert,CEMax=CEMax,Probs=Probs) )

# Free up some memory
rm(Errors,Cert,CEMax,Probs)

# I choose to do this step in R to ensure precision. Highly vectorized
# operations are not very time consuming in R
D[,5:15] <- D[,5:15] / snum
D[,3:4] <- D[,3:4] / cnum

# Add the pattern and parameters to D for clarity
D$pattern <- apply(patmat,2,toString)
D$rmean  <- rm
D$mumean <- um
D$rstd   <- rs
D$mustd  <- us

# Add column labels
colnames(D) <- c("EE","PC","WC","WP","E.0", "E.1", "E.2", "E.3", "E.4", "E.5", "E.6", "E.7", "E.8", "E.9", "E.10","Pattern","rm","rs","um","us")

# Sort by probability, this makes the most sense
D <- D[order(D$PC, decreasing=T),]

# Add row labels
rownames(D) <- 1:nrow(D)

# Get the aggregate means and variances of these statistics
M.EE  <- sum(D$PC*D$EE)
M.WP  <- sum(D$PC*D$WP)
M.WC  <- sum(D$PC*D$WC) 
M.WE0 <- sum(D$PC*D$E.0)

V.EE  <- sum(D$PC*(D$EE - M.EE )^2)
V.WP  <- sum(D$PC*(D$WP - M.WP )^2)
V.WC  <- sum(D$PC*(D$WC - M.WC )^2)
V.WE0 <- sum(D$PC*(D$E.0 - M.WE0)^2)

print(head(D,n=10))

c(  M.EE=M.EE,
	V.EE=V.EE,
	M.WP=M.WP,
	V.WP=V.WP,
	M.WC=M.WC,
	V.WC=V.WC,
	M.WE0=M.WE0,
	V.WE0=V.WE0,
	rm,	rs,	um,	us)

