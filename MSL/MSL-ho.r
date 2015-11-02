# This is the main file to source the data and call the MSL functions

#libraries
library("parallel")
library("pracma")
library("optimx")

cores <- detectCores() - 1

load("../sim/HO.Rda")

#MSL must take in

h <- 100

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

	set.seed(42)

	a <- matrix(c(rnorm(h,mean=rm,sd=rs), rgamma(h,shape=k,scale=t)), ncol=h, byrow=T) 
	b <- data.frame( a )

#	if(anyNA(b[,2])){
#		for(i in 1:1000){
#			b[,2] <- ifelse(is.na(b[,2]) , rgamma(1,shape = k, scale = t),b[,2] )
#			if(!anyNA(b[,2])) break
#		}
#	}
#	if(anyNA(b[,1])){
#		for(i in 1:100){
#			b[,1] <- ifelse(is.na(b[,1]) , rgamma(1,shape = k, scale = t),b[,1] )
#			if(!anyNA(b[,1])) break
#		}
#	}
#	if(any(b[,2]==0)){
#		for(i in 1:1000){
#			b[,2] <- ifelse(b[,2]==0 , rgamma(1,shape = k, scale = t),b[,2] )
#			if(!any(b[,2]==0)) break
#		}
#	}

	sim <- mclapply(b,function(x){

		r  <- 1 - x[1]
		mu <- x[2]

		ctx <- D$max^r/r - D$min^r/r

		UA <- (D$pA0 * D$A0^r/r + D$pA1 * D$A1^r/r + D$pA1 * D$A1^r/r) / ctx / mu
		UB <- (D$pB0 * D$B0^r/r + D$pB1 * D$B1^r/r + D$pB1 * D$B1^r/r) / ctx / mu

		UA.0 <- 0
		UB.0 <- UB - UA

		pA <- exp(UA) / (exp(UA) + exp(UB))
		pB <- exp(UB) / (exp(UA) + exp(UB))

		ll <- ifelse(D$c==0, pA, pB)

		ll

    }, mc.cores=cores)

    prob <- data.frame(matrix(unlist(sim[!is.na(sim)]), ncol=h))
    prob$ID <- D$ID

	llike<-ddply(prob,"ID",function(x){
			x <- x[,!names(x) %in% c("ID")]
				unlist(lapply(x,prod))
			})

	llike<- llike[,!names(llike) %in% c("ID")]

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
#dist <- c(rm=mean(D$r),rs=std(D$r),um=mean(D$mu),us=std(D$mu))

	k <- (mean(D$mu)^2)/(std(D$mu)^2)
	t <- (std(D$mu)^2) / mean(D$mu)

	print(c(k,t))

dist <- c(rm=mean(D$r),rs=std(D$r),k=k,t=t)

dist <- c(rm=mean(D$r),rs=std(D$r),um=mean(D$mu),us=std(D$mu))
init[2:4] <- log(init[2:4])

init <- c(dist)

b <- MSL(init)

lbound <- c(-Inf,0.00000001,.00000001,0.0000001)
ubound <- c(5,1,10,1)

# c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")

#m <- optimx(par=init,fn=MSL, hessian=TRUE)
m <- optim(par=init,fn=MSL, hessian=TRUE, method="BFGS")
#m <- optimx(par=init,fn=MSL, hessian=TRUE, control=(list(all.methods=T)))
#m <- optimx(par=init,fn=MSL, hessian=TRUE, method="L-BFGS-B",lower=lbound,upper=ubound)
#m <- optimx(par=init,fn=MSL, hessian=TRUE, method="L-BFGS-B",lower=lbound,upper=ubound)

#m <- optimx(par=init,fn=MSL, hessian=TRUE, method="CG")
#m <- optim(par=init,fn=MSL, hessian=TRUE, method="SANN")

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


