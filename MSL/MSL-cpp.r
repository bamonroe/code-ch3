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
	h3 <- matrix(H[[3]][,1:h],ncol=h)
	h4 <- matrix(H[[4]][,1:h],ncol=h)

	con <- list(trace=1,maxit=100,REPORT=1)

	# The MSL functions are written entirely in C++ via Rcpp. This cuts down on the
	# optimization time as the Rcpp code is much faster than the native R.

	if(model=="EUT"){
		m <- optim(par=int[1:4],fn=MSL_EUT,method="BFGS", control=con,hessian=TRUE,
				h1=h1,h2=h2,
				A0=D$A0,A1=D$A1,B0=D$B0,B1=D$B1,
				pA0=D$pA0,pA1=D$pA1,pB0=D$pB0,pB1=D$pB1,
				max=D$max,min=D$min,c=D$c
				)
    }
	else if(model=="RDU"){
		m <- optim(par=int[1:8],fn=MSL_RDU,method="BFGS", control=con,hessian=TRUE,
				h1=h1,h2=h2,h3=h3,h4=h4,
				A0=D$A0,A1=D$A1,B0=D$B0,B1=D$B1,
				pA0=D$pA0,pA1=D$pA1,pB0=D$pB0,pB1=D$pB1,
				max=D$max,min=D$min,c=D$c
				)
	}


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
	tt[2:p.num] <- exp(tt[2:p.num])

	ts <- se
	ts[2:p.num] <- exp(se[2:p.num])

	tu <- up
	tu[2:p.num] <- exp(tu[2:p.num])

	tl <- low
	tl[2:p.num] <- exp(tl[2:p.num])

	start <-int[1:p.num]
	start[2:p.num] <- exp(int[2:p.num])

	# Save everything in a convienient place
	mm <- data.frame(real=real[1:p.num],init=start,est=m$par,par=tt,se=se,lower=tl,upper=tu,pvalue=pval,llike=m$value, H=h, HH=HH, UH=UH)
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
dist <- c(.5,.5,.5,.5,1,.5,1,.5)
init <- dist
init[2:8] <- log(dist[2:8])

# do.optim returns a dataset with the relevent information. Let's estimate the
# 4 sets of choice data

load("../data/choice-dat/choice10.Rda")
H <- regenH(D)

real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu),am=1,as=0,bm=1,bs=0)
m10.e <- do.optim(int=init,h=h,model="EUT")
m10.r <- do.optim(int=init,h=h,model="RDU")

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

