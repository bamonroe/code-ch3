# This is the main file to source the data and call the MSL functions

rm(list=ls())

#libraries
library(Rcpp)
library(parallel)
library(microbenchmark)

cores <- detectCores()

# Source the cpp stuff
sourceCpp("../Rcpp/halton.cpp")
sourceCpp("../Rcpp/sim3.cpp")
sourceCpp("../Rcpp/choice.cpp")

# Make a Halton Sequence for use with MSL
#halton(init,H,prime)
mkhalton <- function(h=75,prime=3,HH=2,UH=1,N){

	# 3 ways to determine how many halton draws to make
	H <-ifelse(HH==0,h,ifelse(HH==1,h*N,h*N*T))
	#Halton or RUniform?
	if(UH==1){
		drop <- 0
		hseq <- data.frame(matrix(halton(1+drop,H+drop,prime),ncol=h, byrow=T))
	}
	else{
		hseq <- data.frame(matrix(runif(H),ncol=h, byrow=T))
	}

	# The above will have created either 1 row, N rows, or N*T rows. Need to stack
	# these dataframes on top of themselves until we have N*T rows. But, If we have
	# N rows above, it means we're doing per 'n' draws, so we have to keep them
	# grouped together. So create a number index and sort by this number after the
	# stacking is done.
	hseq$h <- 1:nrow(hseq)

	# How many times do we need to stack?
	stack.num <- nrow(D) / nrow(hseq)

	# Stack 'em up!
	hseq <- matrix( rep( t( hseq ) , stack.num ) , ncol = ncol(hseq) , byrow = TRUE )

	# Now sort by the number index
	hseq<-hseq[order(hseq[,(h+1)]),]

	# Drop the index column
	hseq <- hseq[,1:h]

	return(hseq)
    
}

regenH <- function(Data) {

	# Number of subjects
	N <- max(Data$ID)
	# Number of Tasks
	T <- nrow(Data) / N

	H <- list()

	# Make the sequences
	H[[1]] <- mkhalton(500,prime=3,HH=HH,UH=UH,N=N)
	H[[2]] <- mkhalton(500,prime=7,HH=HH,UH=UH,N=N)

	return(H)
    
}

MSL <- function(par,h1,h2){

	par[2:4] <- exp(par[2:4])

	rm <- par[1]
	rs <- par[2]

	um <- par[3]
	us <- par[4]
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	# Turn the halton sequences into distributions governed by par[]
	r  <- qnorm(h1,mean=rm,sd=rs)
	mu <- qgamma(h2,shape=k,scale=t)

	sim <- simcpp3(r=r,
				mu=mu,
				A0=D$A0,
				A1=D$A1,
				B0=D$B0,
				B1=D$B1,
				pA0=D$pA0,
				pA1=D$pA1,
				pB0=D$pB0,
				pB1=D$pB1,
				max=D$max,
				min=D$min,
				c=D$c 
				)

	# Each column of sim is a simulation, but we need the means of each row,
	# since each row corresponds to a unique choice. Thus the list needs to 
	# be made into a matrix to do sums. 

	#like <- matrix(unlist(sim), ncol=h)

	sl <- sum(log(rowMeans(sim)))

	# optim minimizes functions, so need to return the negitive
	# of the log-likelihood in order to maximize the ll
	return(-sl)

}

do.optim <- function(h,int){
    
	h1 <- matrix(H1[,1:h], ncol=h)
	h2 <- matrix(H2[,1:h], ncol=h)

	con <- list(trace=1,maxit=100)

	m <- optim(par=int,fn=MSL,h=h,h1=h1,h2=h2,method="BFGS", control=con,hessian=TRUE )

	# Get the inverse of the Hessian
	fisher <- solve(m$hessian)
	# Get the square root of it
	se <- sqrt(diag(fisher))
	# Get the 95% confidence interval
	up <- m$par + 1.96*se
	low <- m$par - 1.96*se

	# Get the t stat
	t <- m$par / se

	# Get the p-values
	pval<-2*(1-pt(abs(t),nrow(D)-length(int)))

	# Adjust for the logged parameters
	tt <-  m$par
	tt[2:4] <- exp(tt[2:4])

	ts <- se
	ts[2:4] <- exp(se[2:4])

	tu <- up
	tu[2:4] <- exp(tu[2:4])

	tl <- low
	tl[2:4] <- exp(tl[2:4])

	start <-int
	start[2:4] <- exp(int[2:4])

	# Save everything in a convienient place
	mm <- data.frame(real=real,init=start,est=m$par,par=tt,se=se,lower=tl,upper=tu,pvalue=pval,llike=m$value, H=h, HH=HH, UH=UH)
	# Print these things out
	print(mm)
	cat('\n')
	return(mm)

}


load("Pars.Rda")

## P is a datafram with the population parameters used in the Agg-Gen
#load("../data/agg-dat/P-Agg30k-S500.Rda")
#
#
## Now each column is a set
#P <- data.frame(t(data.matrix(P)))
#
#DD <- rbind( P[1:2,] , (P[3,]^2 / P[4,]^2) , P[4,]^2 / P[3,] )
#
#DD <- DD[,which(DD[1,] > -1.8 & DD[1,] < 1.3)]
#
#save(DD,file="Pars.Rda")
#
#rm(P)



# The instrument
A0 <- 1.60
A1 <- 2.00
B0 <- .10
B1 <- 3.85

pA0 <- seq(from=0.9, to=0, by= -.1)
pB0 <- seq(from=0.9, to=0, by= -.1)
pA1 <- seq(from=0.1, to=1, by= .1)
pB1 <- seq(from=0.1, to=1, by= .1)

A <- matrix(c(A0,A1),nrow=10,ncol=2,byrow=T)
B <- matrix(c(B0,B1),nrow=10,ncol=2,byrow=T)

pA <- matrix(c(pA0,pA1),ncol=2,nrow=10,byrow=F)
pB <- matrix(c(pB0,pB1),ncol=2,nrow=10,byrow=F)

Max <- rep(3.85,10)
Min <- c(rep(0.1,9),2)

# Draw 's.num' samples from the distribution, have each sample complete the
# task 't.num' times 


# Make a matrix of R values and Mu values from the populations

subjects <- 150
s.num    <- 10

A   <- do.call(rbind, replicate(subjects,A,simplify=F))
B   <- do.call(rbind, replicate(subjects,B,simplify=F))
pA  <- do.call(rbind, replicate(subjects,pA,simplify=F))
pB  <- do.call(rbind, replicate(subjects,pB,simplify=F))
Max <- do.call(rbind, replicate(subjects,Max,simplify=F))
Min <- do.call(rbind, replicate(subjects,Min,simplify=F))


r.mat  <- matrix( rnorm( subjects* s.num ,mean=DD[1,7],sd=DD[2,7]), ncol=s.num )
mu.mat <- matrix( rgamma(subjects* s.num ,shape=DD[3,7],scale=DD[4,7]),ncol=s.num )

r.mat  <- matrix( rnorm( subjects* s.num ,mean=.45,sd=.1), ncol=s.num )
mu.mat <- matrix( rgamma(subjects* s.num ,shape=14,scale=.01),ncol=s.num )

r.mat  <- cbind( 1:nrow(r.mat) , r.mat)
mu.mat <- cbind( 1:nrow(mu.mat) , mu.mat)

r.mat  <- do.call(rbind, replicate(10,r.mat,simplify=F))
mu.mat <- do.call(rbind, replicate(10,mu.mat,simplify=F))

r.mat  <- r.mat[order(r.mat[,1]),2:ncol(r.mat)]
mu.mat <- mu.mat[order(mu.mat[,1]),2:ncol(mu.mat)]

par.mat <- data.frame(rbind(r.mat,mu.mat))

# getChoice returns a matrix of choices given the instrument and the
# preferences

t.num <- 100

mkChoice <- function(par.vec){

	r.vec <- par.vec[1:(length(par.vec)/2)]
	mu.vec <- par.vec[(length(par.vec)/2)+1:length(par.vec)]

	getChoice(r.vec,mu.vec,A,B,pA,pB,Max,Min,t.num)
}


OO <- mclapply(par.mat,mkChoice,mc.cores=cores)
	

stop("here")


# How many simulations to be done?
h <- 100

# Which type of halton draws? 0 means every observation gets the same draw, 1
# means that each subject gets their own draw, 2 means every observaltion gets
# different draw.
HH <- 2

# Use the halton sequences (1) or use the urand() function (2) to get uniform
# distributions?
UH <- 1

# What are the initial values optim will start with? Log the ones that need to
# be greater than 0.
dist <- c(.5,.5,.5,.5)
init <- dist
init[2:4] <- log(dist[2:4])

# do.optim returns a data.frame with the relevent information. 

load("../data/choice-dat/choice10.Rda")
H <- regenH(D)
H1 <- H[[1]]
H2 <- H[[2]]
real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))


	m10 <- do.optim0(int=init,h=h )


