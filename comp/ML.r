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
	mu <- par[2]

	UA <- dd$pA0 * rr(dd$A0,r) + dd$pA1 * rr(dd$A1,r)
	UB <- dd$pB0 * rr(dd$B0,r) + dd$pB1 * rr(dd$B1,r)

	ctx <- rr(dd$max,r) - rr(dd$min,r)

	UB <- UB/ctx/mu - UA/ctx/mu
	UA <- 0
	
	eA <- exp(UA)
	eB <- exp(UB)
	ee <- eA + eB

	pA <- eA / ee
	pB <- eB / ee

	ll <- ifelse(D$c==0, log(pA), log(pB))

	-sum(ll)
}

ML.RDU <- function(par){

	dd <- E

	r <- par[1]
	mu <- par[2]
	a <- par[3]
	b <- par[4]

	wA1 <- pw(dd$pA1,a,b)
	wA0 <- pw(dd$pA0 + dd$pA1,a,b) - pw(dd$pA1,a,b)

	wB1 <- pw(dd$pB1,a,b)
	wB0 <- pw(dd$pB0 + dd$pB1,a,b) - pw(dd$pB1,a,b)

	UA <- wA0 * rr(dd$A0,r) + wA1 * rr(dd$A1,r) 
	UB <- wB0 * rr(dd$B0,r) + wB1 * rr(dd$B1,r) 

	ctx <- rr(dd$max,r) - rr(dd$min,r)

	UB <- UB/ctx/mu 
	US <- UA/ctx/mu

	UB <- UB - UA
	UA <- 0

	eA <- exp(UA)
	eB <- exp(UB)
	ee <- eA + eB

	pA <- eA / ee
	pB <- eB / ee

	ll <- ifelse(dd$c==0, log(pA), log(pB))

	-sum(ll)
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
	i1 <- grepl('^p', names(tmp))

	tmp[] <- t(apply(tmp, 1, function(x) {i2 <- order(x[i1]*x[!i1])
								c(x[i1][i2], x[!i1][i2])}))
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
init <- c(r=.65,mu=.3,a=1,b=1)
RDU  <- optim(par=init,fn=ML.RDU,hessian=TRUE)

init <- c(r=.65,mu=.3)
EUT  <- optim(par=init,fn=ML.EUT,hessian=TRUE)

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

# Print these things out
print(RDU$par)
print(R.se)
print(R.up)
print(R.low)

print(EUT$par)
print(E.se)
print(E.up)
print(E.low)


