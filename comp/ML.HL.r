#load("choice.Rda")

load("choice.Rda")

# optim minimizes functions, so need to return the negitive
# of the log-likelihood in order to maximize the ll

pw <- function(p,a,b){
    exp( -b * (-log(p))^a )
}

rr <- function(x,r){
    
    ifelse(r == 1 , log(x),	x^(1-r) / (1-r) )
}


ML.EUT <- function(par){

	dd <- E

	r <- par[1]
	mu <- exp( par[2])

	UA <- dd$pA0 * rr(dd$A0,r) + dd$pA1 * rr(dd$A1,r)
	UB <- dd$pB0 * rr(dd$B0,r) + dd$pB1 * rr(dd$B1,r)

	ctx <- sapply(dd$max,rr,r=r) - sapply(dd$min,rr,r=r)

	UB <- UB/ctx/mu 
	UA <- UA/ctx/mu

	UB <- UB - UA
	UA <- 0

	eA <- exp(UA)
	eB <- exp(UB)
	ee <- eA + eB

	pA <- eA / ee
	pB <- eB / ee

	ll <- ifelse(D$c==0, log(pA), log(pB))

	-sum(ll,na.rm=T)

}

ML.RDU <- function(par){

	dd <- E

	r <- par[1]
	mu <- exp( par[2])
	a <- par[3]
	b <- par[4]

	wA1 <- pw(dd$pA1,a,b)
	wA0 <- pw(dd$pA0 + dd$pA1,a,b) - pw(dd$pA1,a,b)

	wB1 <- pw(dd$pB1,a,b)
	wB0 <- pw(dd$pB0 + dd$pB1,a,b) - pw(dd$pB1,a,b)

	UA <- wA0 * rr(dd$A0,r) + wA1 * rr(dd$A1,r) 
	UB <- wB0 * rr(dd$B0,r) + wB1 * rr(dd$B1,r) 

	ctx <- sapply(dd$max,rr,r=r) - sapply(dd$min,rr,r=r)

	UB <- UB/ctx/mu 
	UA <- UA/ctx/mu

	UB <- UB - UA
	UA <- 0

	eA <- exp(UA)
	eB <- exp(UB)
	ee <- eA + eB

	pA <- eA / ee
	pB <- eB / ee

	ll <- ifelse(dd$c==0, log(pA), log(pB))

	-sum(ll,na.rm=T)
}

ML.MIX <- function(par){

	dd <- E

	pi <- par[7]

	mEUT <- exp(pi) / ( exp(pi) + 1 )
	mRDU <- 1 - mEUT

	r <- par[1]
	mu <- exp( par[2])
	a <- par[3]
	b <- par[4]

	wA1 <- pw(dd$pA1,a,b)
	wA0 <- pw(dd$pA0 + dd$pA1,a,b) - pw(dd$pA1,a,b)

	wB1 <- pw(dd$pB1,a,b)
	wB0 <- pw(dd$pB0 + dd$pB1,a,b) - pw(dd$pB1,a,b)

	UA <- wA0 * rr(dd$A0,r) + wA1 * rr(dd$A1,r) 
	UB <- wB0 * rr(dd$B0,r) + wB1 * rr(dd$B1,r) 

	ctx <- sapply(dd$max,rr,r=r) - sapply(dd$min,rr,r=r)

	UB <- UB/ctx/mu 
	UA <- UA/ctx/mu

	UB <- UB - UA
	UA <- 0

	eA <- exp(UA)
	eB <- exp(UB)
	ee <- eA + eB

	pA <- eA / ee
	pB <- eB / ee

	pRDU <- ifelse(dd$c==0, pA, pB)


	# Now for EUT

	r <-  par[5]
	mu <- exp( par[6])

	UA <- dd$pA0 * rr(dd$A0,r) + dd$pA1 * rr(dd$A1,r)
	UB <- dd$pB0 * rr(dd$B0,r) + dd$pB1 * rr(dd$B1,r)

	ctx <- sapply(dd$max,rr,r=r) - sapply(dd$min,rr,r=r)

	UB <- UB/ctx/mu 
	UA <- UA/ctx/mu

	UB <- UB - UA
	UA <- 0

	eA <- exp(UA)
	eB <- exp(UB)
	ee <- eA + eB

	pA <- eA / ee
	pB <- eB / ee

	pEUT <- ifelse(D$c==0, pA, pB)

	# The nuputal

	j <- log( mEUT * pEUT + mRDU * pRDU )


	ll<-(-sum(j,na.rm=T))

	print(ll)

	ll

    
}

# Rank non-zero probability outcomes and their probabilities
# Credit to user "akrun" from www.stackoverflow.com for providing a very simple
# method to do the ranking. This saved me a fair amount of time and headache
# and makes the code far more readable and flexible than the solution I was
# planning on using.

# Question and answer:
#http://stackoverflow.com/questions/32598517/how-to-rearrange-elements-in-rows-of-a-dataframe-in-r-based-on-conditions


rank <- function(i){
    
	gr <- paste("[",i,"]\\d{1}",sep="")
	tmp <-  D[grepl(gr,names(D))]

	# The number of columns divided by 2 is the number of outcomes
	n <- ncol(tmp) / 2

	tmp[] <- t(apply(tmp, 1, function(x) {

				# x is the row , first n elements are probs, the second n
				# elements are the corresponding outcomes

				uo <- c()	# vector for unordered outcomes
				up <- c()   # vector for unordered probabilities
				oo <- c()   # vector for ordered outcomes
				op <- c()   # vector for ordered probabilities

				for (i in 1:n){				# Loop through probabilities
					if( x[i] != 0){			# if probability isn't 0, it needs to be ordered
                        op <- c(op, x[i])	# add the probability to the list
                        oo <- c(oo, x[i+n]) # add the outcome to the list
					}
					else{					# if the probability is 0, add it isn't ordered
                        up <- c(up, x[i] )  
                        uo <- c(uo, x[i+n] )
					}
				}

				r <- order(oo)	# Order the elements of the outcomes list

				p <- c(up, op[r]) # the vector of probabilites with the 0's at the back
				o <- c(uo, oo[r]) # vector of outcomes with 0 probability outcomes in the back

				c(p,o)

			}))

	return(tmp)

}

opt <- c("A","B")
E <- rank(opt[1])

if (length(opt) > 1) {
	for (i in 2:length(opt)){
		E <- cbind(E, rank(opt[i]))
	}
}

E$ID <- D$ID
E$c <- D$c
E$max <- D$max
E$min <- D$min

# To get the raw output, just run like this
# you can see that $par contains the results,
# value contains the LL, etc.

con <- list(trace=2)
con <- list()

rinit <- c(r=.65,mu=log(.3),a=1,b=1)

ML.RDU(rinit)

#stop("here")

#RDU  <- optim(par=rinit,fn=ML.RDU,hessian=TRUE,control=con)

einit <- c(r=.65,mu=log(.3))

#EUT  <- optim(par=einit,fn=ML.EUT,hessian=TRUE,control=con)

minit <- c(RDU=rinit,EUT=einit,pi=-1)

MIX  <- optim(par=minit,fn=ML.MIX,hessian=TRUE,control=con)

# We need to get the standard errors of the parameters.
# The standard error is the square root of the the inverse
# of the identity of the negitive hessian matrix. Since we
# were minimizing the negitive of the log-likilhood, we
# are returned the negitive hessian

# Get the inverse of the Hessian
R.fisher <- solve(RDU$hessian)
# Get the square root of it
R.se <- sqrt(diag(R.fisher))
# Get the 95% confidence interval
R.up <- RDU$par + 1.96*R.se
R.low <- RDU$par - 1.96*R.se

# Get the inverse of the Hessian
E.fisher <- solve(EUT$hessian)
# Get the square root of it
E.se <- sqrt(diag(E.fisher))
# Get the 95% confidence interval
E.up  <- EUT$par + 1.96*E.se
E.low <- EUT$par - 1.96*E.se

# Get the inverse of the Hessian
M.fisher <- solve(MIX$hessian)
# Get the square root of it
M.se <- sqrt(diag(M.fisher))
# Get the 95% confidence interval
M.up  <- MIX$par + 1.96*M.se
M.low <- MIX$par - 1.96*M.se

# Print these things out
print(RDU$par)
print(R.se)
print(R.up)
print(R.low)

print(EUT$par)
print(E.se)
print(E.up)
print(E.low)

print(MIX$par)
print(M.se)
print(M.up)
print(M.low)

