# This is the main file to source the data and call the MSL functions

#libraries
library("parallel")
library("plyr")

cores <- detectCores() -1

load("../sim/choice.Rda")

# Number of subjects
N <- max(D$ID)

# Number of Tasks
T <- nrow(D) / N

# Make a Halton Sequence for use with MSL
halton <- function(init,H,prime){
	sapply(c(init:H),function(index,base=prime){
    
		r <- 0
		f <- 1
		i <- index

		while(i > 0){
			f <- f / base
			r <- r + f*(i %% base)
			i <- floor(i /base)
		}
		return(r)
	})
}

logit.normal <- function(x,m,s) {
    a <- 1 / (s * sqrt(2*pi) )
    b <- 1 / (x * (1-x))
    c <- exp( -( (log(x/(1-x)) - m)^2 / (2 * s^2) ) )
	a*b*c
}

# How many simulations will be performed
# With B^h scenerio, this is the same as the size of the sequnce needed
h <- 150

# With B^h_n scenerio, h needs to be multiplied by N to get the sequence
# size
NH <- h*N

# With B^h_ni scenerio, h needs to be multiplied by N and T to get the sequence
# size
NHT <- h*N*T

# Need to make a different halton sequence for each distribution
H <- NH
# Toss the first few values of the sequence because they're ever so slightly
# correlated
drop <- 30

h1 <- data.frame(matrix(halton(1+drop,H+drop,3),ncol=h))
h2 <- data.frame(matrix(halton(1+drop,H+drop,7),ncol=h))

# The above will have created either 1 row, N rows, or N*T rows. Need to stack
# these dataframes on top of themselves until we have N*T rows. But, If we have
# N rows above, it means we're doing per 'n' draws, so we have to keep them
# grouped together. So create a number index and sort by this number after the
# stacking is done.
h1$h <- 1:nrow(h1)
h2$h <- 1:nrow(h2)

# How many times do we need to stack?
stack.num <- nrow(D) / nrow(h1)

# Stack 'em up!
h1 <- matrix( rep( t( h1 ) , stack.num ) , ncol = ncol(h1) , byrow = TRUE )
h2 <- matrix( rep( t( h2 ) , stack.num ) , ncol = ncol(h2) , byrow = TRUE )

# Now sort by the number index
h1<-hr[order(h1[,(h+1)]),]
h2<-hm[order(h2[,(h+1)]),]

# Drop the index column
h1 <- h1[,1:h]
h2 <- h2[,1:h]

stop("here")

MSL <- function(par){

	par[2:4] <- exp(par[2:4])

	rm <- par[1]
	rs <- par[2]

	um <- par[3]
	us <- par[4]
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	#k <- par[3]
	#t <- par[4]

	# Turn the halton sequences into distributions governed by par[]
	a <- matrix(qnorm(h1,mean=rm,sd=rs))
	b <- matrix(qgamma(h2,shape=k,scale=t)) 

	dd <- data.frame( rbind(a,b) )

	sim <- mclapply(dd,function(x){

		r  <- 1 - x[1:N]
		mu <- x[(N+1):(2*N)]

		ctx <- D$max^r/r - D$min^r/r

		UA <- (D$pA0 * D$A0^r/r + D$pA1 * D$A1^r/r) 
		UB <- (D$pB0 * D$B0^r/r + D$pB1 * D$B1^r/r) 

		# Calculte the probability stuff, need to rebase the utility because 
		# of precision and machine limit problems. exp(>709) evalutes as Inf.
		
		UB.1 <- UB/ctx/mu - UA/ctx/mu
		UA.1 <- 0

		# Have things gone haywire because of the computer's inability to handle
		# numbers bigger than ~3e310 ?
		c.N <- is.nan(UB.1) | is.infinite(UB.1)

		pA <- ifelse( c.N ,		# Are we dealing with an insane number?
			# yes:
			ifelse( UB > UA , 0 , 1 ) ,
			# no, but are we making an insane number via exp?
			ifelse( UB.1 > 709 , 0 , 
				ifelse( UB.1 < -709, 1 , exp(UA.1) / (exp(UA.1) + exp(UB.1)) )
			)
		)

		pB <- 1 - pA

		ll<-ifelse(D$c==0,pA,pB)
		ll

    }, mc.cores=cores)

    prob <- data.frame(matrix(unlist(sim[!is.na(sim)]), ncol=h))
	llike <-prob
#    prob$ID <- D$ID

#    return(prob)

#	llike<-ddply(prob,"ID",function(x){
#			x <- x[,!names(x) %in% c("ID")]
#				unlist(lapply(x,prod))
#			})

#	llike<- llike[,!names(llike) %in% c("ID")]

	sl <- rowMeans(llike, na.rm=T)
	sl <- sum(log(sl))
	print(par)
	print(sl)
	return(-sl)

}

# optim minimizes functions, so need to return the negitive
# of the log-likelihood in order to maximize the ll

# To get the raw output, just run like this
# you can see that $par contains the results,
# value contains the LL, etc.

# Assign the results to an object for recalling later
#dist <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))

	k <- (mean(D$mu)^2)/(sd(D$mu)^2)
	t <- (sd(D$mu)^2) / mean(D$mu)

	print(c(k,t))

dist <- c(rm=mean(D$r),rs=sd(D$r),k=k,t=t)

dist <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))
#dist <- c(rm=rnorm(1),rs=runif(1),um=runif(1),us=runif(1))

init <- c(dist)

init <- c(rm=0.5984992,rs=0.1047776,um=0.3572698,us= 0.1025885 )

init[2:4] <- log(init[2:4])
test <- MSL(init)

#stop("this")

lbound <- c(-Inf,0.00000001,.00000001,0.0000001)
ubound <- c(Inf,Inf,Inf,Inf)

# c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")

m <- optim(par=init,fn=MSL, hessian=TRUE, method="Nelder-Mead")
#m <- optim(par=init,fn=MSL, hessian=TRUE, method="BFGS")

# We need to get the standard errors of the parameters.
# The standard error is the square root of the the inverse
# of the identity of the negitive hessian matrix. Since we
# were minimizing the negitive of the log-likilhood, we
# are returned the negitive hessian

# Get the inverse of the Hessian
fisher <- solve(m$hessian)
# Get the square root of it
se <- sqrt(diag(fisher))
# Get the 95% confidence interval
up <- m$par + 1.96*se
low <- m$par - 1.96*se

# Print these things out
print(fisher)
print(m$par)
print(se)
print(up)
print(low)


