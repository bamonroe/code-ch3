# optim minimizes functions, so need to return the negitive
# of the log-likelihood in order to maximize the ll

ML.EUT <- function(par){

	r <- par[1]
	mu <- par[2]

	UA <- D$pA0 * CRRA(D$A0,r) + D$pA1 * CRRA(D$A1,r)
	UB <- D$pB0 * CRRA(D$B0,r) + D$pB1 * CRRA(D$B1,r)

	pA <- exp(UA/mu) / (exp(UA/mu) + exp(UB/mu))
	pB <- exp(UB/mu) / (exp(UA/mu) + exp(UB/mu))

	ll <- ifelse(D$c==0, log(pA), log(pB))

	-sum(ll)
}

# To get the raw output, just run like this
# you can see that $par contains the results,
# value contains the LL, etc.
optim(c(r=.65,mu=.3),fn=ML.EUT, hessian=TRUE)

# Assign the results to an object for recalling later
m <- optim(c(r=.65,mu=.3),fn=ML.EUT, hessian=TRUE)

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
print(m$par)
print(se)
print(up)
print(low)

