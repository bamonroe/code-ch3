# This is the main file to source the data and call the MSL functions

#libraries

library("microbenchmark")
library("rbenchmark")

library("parallel")
cores <- detectCores()

library(Rcpp)
# Source the cpp stuff
sourceCpp("~/git/thesis/halton/halton.cpp")
sourceCpp("~/git/thesis/Rcpp/simply.cpp")

# I use openmp in the cpp  files so I need to add additional compiler flags
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

sourceCpp("~/git/thesis/Rcpp/sim2.cpp")
sourceCpp("~/git/thesis/Rcpp/sim3.cpp")
sourceCpp("~/git/thesis/Rcpp/sim4.cpp")

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

MSL1 <- function(par,h,h1,h2){

	par[2:4] <- exp(par[2:4])

	rm <- par[1]
	rs <- par[2]

	um <- par[3]
	us <- par[4]
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	# Turn the halton sequences into distributions governed by par[]
	a <- matrix(qnorm(h1,mean=rm,sd=rs), ncol=h)
	b <- matrix(qgamma(h2,shape=k,scale=t), ncol=h) 

	dd <- data.frame( rbind(a,b) )

	sim <- mclapply(dd, sim1 ,
	#sim <- lapply(dd, simcpply ,
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
						c=D$c, 
						#)
						mc.cores=cores)

	# Each column of sim is a simulation, but we need the means of each row,
	# since each row corresponds to a unique choice. Thus the list needs to 
	# be made into a matrix to do sums. 

	like <- matrix(unlist(sim), ncol=h)
	sl <- sum(log(rowMeans(like)))

	# optim minimizes functions, so need to return the negitive
	# of the log-likelihood in order to maximize the ll
	return(-sl)

}

MSL2 <- function(par,h,h1,h2){

	par[2:4] <- exp(par[2:4])

	rm <- par[1]
	rs <- par[2]

	um <- par[3]
	us <- par[4]
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	# Turn the halton sequences into distributions governed by par[]
	a <- matrix(qnorm(h1,mean=rm,sd=rs), ncol=h)
	b <- matrix(qgamma(h2,shape=k,scale=t), ncol=h) 

	dd <- rbind(a,b) 

	sim <- sim2(DD=dd, 
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

MSL3 <- function(par,h,h1,h2){

	par[2:4] <- exp(par[2:4])

	rm <- par[1]
	rs <- par[2]

	um <- par[3]
	us <- par[4]
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	# Turn the halton sequences into distributions governed by par[]
	a <- matrix(qnorm(h1,mean=rm,sd=rs), ncol=h)
	b <- matrix(qgamma(h2,shape=k,scale=t), ncol=h) 

	dd <- rbind(a,b) 

	sim <- sim3(DD=dd, 
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

MSL4 <- function(par,h,h1,h2){

	par[2:4] <- exp(par[2:4])

	rm <- par[1]
	rs <- par[2]

	um <- par[3]
	us <- par[4]
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	# Turn the halton sequences into distributions governed by par[]
	a <- matrix(qnorm(h1,mean=rm,sd=rs), ncol=h)
	b <- matrix(qgamma(h2,shape=k,scale=t), ncol=h) 

	dd <- rbind(a,b) 

	sim <- sim4(DD=dd, 
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

do.optim <- function(h,int,Model){
    
	h1 <- H1[,1:h]
	h2 <- H2[,1:h]

	con <- list(trace=1,maxit=100)

	m <- optim(par=int,fn=Model,h=h,h1=h1,h2=h2,method="BFGS", control=con,hessian=TRUE )

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
	print("")
	return(mm)

}

# How many simulations to be done?
h <- 250

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

# do.optim returns a dataset with the relevent information. Let's estimate the
# 4 sets of choice data

load("../sim/choice10.Rda")
H <- regenH(D)
H1 <- H[[1]]
H2 <- H[[2]]
real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))

h1 <- H1[,1:h]
h2 <- H2[,1:h]


#for(i in 1:5){
#    
#
#s1 <- proc.time()
#mm0 <- MSL0(par=init,h=h,h1=h1,h2=h2)
#e1 <- proc.time()
#
#s2 <- proc.time()
#mm2 <- MSL2(par=init,h=h,h1=h1,h2=h2)
#e2 <- proc.time()
#
#t1 <- e1-s1
#t2 <- e2-s2
#
#print(t1-t2)
#
#}

mm1 <- MSL1(par=init,h=h,h1=h1,h2=h2)
mm2 <- MSL2(par=init,h=h,h1=h1,h2=h2)
mm3 <- MSL3(par=init,h=h,h1=h1,h2=h2)
mm4 <- MSL4(par=init,h=h,h1=h1,h2=h2)



print(identical(mm1,mm2,mm3,mm4))

print(c(mm1,mm2,mm3,mm4))


#microbenchmark(
#	mm1 <- MSL1(par=init,h=h,h1=h1,h2=h2),
#	mm2 <- MSL2(par=init,h=h,h1=h1,h2=h2),
#	mm3 <- MSL3(par=init,h=h,h1=h1,h2=h2),
#	mm4 <- MSL4(par=init,h=h,h1=h1,h2=h2),
#	times=100
#)

benchmark(
	mm1 <- MSL1(par=init,h=h,h1=h1,h2=h2),
	mm2 <- MSL2(par=init,h=h,h1=h1,h2=h2),
	mm3 <- MSL3(par=init,h=h,h1=h1,h2=h2),
	mm4 <- MSL4(par=init,h=h,h1=h1,h2=h2),
	replications=10
)


stop("here")

#load("../sim/choice20.Rda")
#H <- regenH(D)
#H1 <- H[[1]]
#H2 <- H[[2]]
#real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))
#m20 <- do.optim(int=init,h=h )
#
#load("../sim/choice30.Rda")
#H <- regenH(D)
#H1 <- H[[1]]
#H2 <- H[[2]]
#real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))
#m30 <- do.optim(int=init,h=h )
#
#load("../sim/choice40.Rda")
#H <- regenH(D)
#H1 <- H[[1]]
#H2 <- H[[2]]
#real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))
#m40 <- do.optim(int=init,h=h )



