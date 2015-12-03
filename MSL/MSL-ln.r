# This is the main file to source the data and call the MSL functions

#libraries
library(Rcpp)

# Source the cpp stuff
sourceCpp("../Rcpp/halton.cpp")
sourceCpp("../Rcpp/sim.cpp")

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
	H[[3]] <- mkhalton(500,prime=13,HH=HH,UH=UH,N=N)
	H[[4]] <- mkhalton(500,prime=17,HH=HH,UH=UH,N=N)

	return(H)
    
}

do.optim <- function(int,h=50,model="EUT"){
    
	h1 <- matrix(H[[1]][,1:h],ncol=h)
	h2 <- matrix(H[[2]][,1:h],ncol=h)

	con <- list(trace=1,maxit=1000,REPORT=1)

	# The MSL functions are written entirely in C++ via Rcpp. This cuts down on the
	# optimization time as the Rcpp code is much faster than the native R.

	if(model=="EUT"){
		m <- optim(par=int[1:7],fn=MSL_EUT_LN,method="BFGS", control=con,hessian=TRUE,
				h1=h1,h2=h2,
				A=cbind(D$A0,D$A1),B=cbind(D$B0,D$B1),
				pA=cbind(D$pA0,D$pA1),pB=cbind(D$pB0,D$pB1),
				max=D$max,min=D$min,c=D$c
				)
    }
	else
		return(1)


	# p.num is the number of parameters
	p.num <- length(m$par)

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
	pval<-2*(1-pt(abs(t),nrow(D)-p.num))

	# Adjust for the logged parameters
	tt <-  m$par
	tt <- c(tt[1],exp(tt[2:3]),tt[4:5],exp(tt[6:7]))

	ts <- se
	ts <- c(ts[1],exp(ts[2:3]),ts[4:5],exp(ts[6:7]))

	tu <- up
	ts <- c(ts[1],exp(ts[2:3]),ts[4:5],exp(ts[6:7]))

	tl <- low
	tl <- c(tl[1],exp(tl[2:3]),tl[4:5],exp(tl[6:7]))

	start <-int[1:p.num]
	start <- c(start[1],exp(start[2:3]),start[4:5],exp(start[6:7]))

	# Save everything in a convienient place
	if(exists("real")){
        
		mm <- data.frame(real=real[1:p.num],init=start,est=m$par,par=tt,se=se,lower=tl,upper=tu,pvalue=pval,llike=m$value, H=h, HH=HH, UH=UH)
	}
	else{
		mm <- data.frame(init=start,est=m$par,par=tt,se=se,lower=tl,upper=tu,pvalue=pval,llike=m$value, H=h, HH=HH, UH=UH)
	}

	# Print these things out
	print(mm)
	cat('\n')
	return(mm)

}


# How many simulations to be done?
h <- 300

# Which type of halton draws? 0 means every observation gets the same draw, 1
# means that each subject gets their own draw, 2 means every observaltion gets
# different draw.
HH <- 2

# Use the halton sequences (1) or use the urand() function (2) to get uniform
# distributions?
UH <- 1

# What are the initial values optim will start with? Log the ones that need to
# be greater than 0.
c(0.007813979,  0.177822618,  7.573731819, -3.190394042)
dist <- c(-2,.5,15,0,-4,.6,14)
init <- c(dist[1],log(dist[2]),log(dist[3]),dist[4],dist[5],log(dist[6]),log(dist[7]))

# do.optim returns a dataset with the relevent information. Let's estimate the
# 4 sets of choice data

load("../data/choice-dat/choice10.Rda")
H <- regenH(D)

real0 <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))
m10.e <- do.optim(int=init,h=h,model="EUT")
#m10.r <- do.optim(int=init,h=h,model="RDU")

#load("../sim/choice20.Rda")
#H <- regenH(D)
#real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))
#m20 <- do.optim(int=init,h=h )
#
#load("../sim/choice30.Rda")
#H <- regenH(D)
#real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))
#m30 <- do.optim(int=init,h=h )
#
#load("../sim/choice40.Rda")
#H <- regenH(D)
#real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))
#m40 <- do.optim(int=init,h=h )

