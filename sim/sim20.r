# This whole bit works, but I'm certain its not the most efficient way to
# do this.

#Using the Holt Laury(2002) instrument
labels <- list("A0","A1","B0","B1","pA0","pA1","pB0","pB1")
numbers <- c(1.60, 2.00, 0.10, 3.85, .9, .1, .9, .1)

subjects <- 150
N <- subjects * 10

A0 <- 1.60
A1 <- 2.00
B0 <- .10
B1 <- 3.85

pA0 <- seq(from=0.9, to=0, by= -.1)
pB0 <- seq(from=0.9, to=0, by= -.1)
pA1 <- seq(from=0.1, to=1, by= .1)
pB1 <- seq(from=0.1, to=1, by= .1)

# Define some population parameters, mean and standard devitions of r
# and fechner values
rm <- 0.65
rs <- 0.35

um <- 0.35
us <- 0.3

# Fechner will use a gamma distribution, so need to back out shape and
# scale parameters
k <- (um^2)/(us^2)
t <- (us^2)/um

# grab one value from this distribution
rval <-rnorm(1,mean = rm, sd = rs)
uval <-rgamma(1,shape = k, scale = t)

# R can have multiple datasets open at once. It calls these "data frames"
# I can think of a lot of instances where this would be hugely beneficial
d <- data.frame(pA0=pA0,pA1=pA1,pB0=pB0,pB1=pB1,A0=A0,A1=A1,B0=B0,B1=B1)

rm(A0,A1,B0,B1,pA0,pA1,pB0,pB1)

# For Use with context
d$min <- ifelse(d$pB0 > 0,d$B0,d$A1)
d$max <- d$B1

# rbind d multiple times to increase the number of MPLs
d <- rbind(d,d)

d$lotid <- 1:nrow(d)

D <- d
D$ID <- 1
D$r  <-rnorm(1,mean = rm, sd = rs)
D$mu <-rgamma(1,shape = k, scale = t) 


for(i in 2:subjects){

	d$ID <- i
	d$r  <-rnorm(1,mean = rm, sd = rs)
	d$mu <-rgamma(1,shape = k, scale = t)

	D <- rbind(D,d)
}

# Practice using R as a functional language. This is a cheap example
# But you cannot do the equivaluent of this in Stata with "programs"
CRRA <- function(out,r){
	((out^(1-r)) / (1-r))
}

# Adding the Utility of A and B to the dataframe, this isn't actually
# necessary as it is in R. We could have just defined UA and UB as 
# vectors and they would have assumed the correct length. The benefit of
# not putting them in the dataframe would be that when the dataframe is
# saved, we don't keep these intermediary variables

ctx <- CRRA(D$max,D$r) - CRRA(D$min,D$r)

D$UA <- (D$pA0 * CRRA(D$A0,D$r) + D$pA1 * CRRA(D$A1,D$r))
D$UB <- (D$pB0 * CRRA(D$B0,D$r) + D$pB1 * CRRA(D$B1,D$r))

UB.1 <- (D$UB/ctx/D$mu) - (D$UA/ctx/D$mu)
UA.1 <- 0

# Have things gone haywire because of the computer's inability to handle
# numbers bigger than ~3e310 ?
c.N <- is.nan(UB.1) | is.infinite(UB.1)

D$pA <- ifelse( c.N ,		# Are we dealing with an insane number?
	# yes:
	ifelse( D$UB > D$UA , 0 , 1 ) ,
	# no, but are we making an insane number via exp?
	ifelse( UB.1 > 709 , 0 , 
		ifelse( UB.1 < -709, 1 , exp(UA.1) / (exp(UA.1) + exp(UB.1)) )
	)
)

D$pB <- 1 - D$pA

# Random uniform number
rand <- runif(nrow(D))

# This is a great R function, ifelse collapses would would otherwise
# potentially be several lines of code into a very readable one line
# statement.
D$c <- ifelse(D$pA > rand, 0, 1)

save(D ,file="choice20.Rda")

# Dataset is generated, call the file that does ML estimation.
#source("ML.r", echo=TRUE,print.eval=TRUE)
