# This whole bit works, but I'm certain its not the most efficient way to
# do this.

#Using the Holt Laury(2002) instrument
labels <- list("A0","A1","B0","B1","pA0","pA1","pB0","pB1")
numbers <- c(1.60, 2.00, 0.10, 3.85, .9, .1, .9, .1)

subjects <- 150

A0 <- c(rep(1,10),rep(10,15))
A1 <- c(rep(20,10),rep(34,15))
A2 <- c(rep(41,10),rep(51,15))
B0 <- c(rep(10,10),rep(20,15))
B1 <- c(rep(20,10),rep(34,15))
B2 <- c(rep(30,10),rep(41,15))

pA0 <- c(0.625,0.375,0,0.125,0.5,0.25,0.25,0.25,0.125,0.125,0.125,0.25,0.625,0.125,0.125,0.375,0,0.5,0.75,0.25,0,0,0.25,0.5,0.25 )
pA1 <- c(0,0.625,1,0.75,0.375,0.75,0.625,0.25,0.375,0.25,0.875,0.75,0.375,0.5,0.75,0.375,0.75,0.125,0,0.375,0.875,0.625,0.5,0.5,0.5	) 
pA2 <- c(.375,0,0,0.125,0.125,0,0.125,0.5,0.5,0.625,0,0,0,0.375,0.125,0.25,0.25,0.375,0.25,0.375,0.125,0.375,0.25,0,0.25 )

pB0 <- c(0.375,0.5,0.125,0.25,0.625,0.375,0.375,0.125,0,0,0.25,0.5,0.75,0.25,0.375,0.5,0.125,0.375,0.625,0.375,0.125,0.125,0.125,0.625,0.375 )
pB1 <- c(0.625,0.25,0.5,0.5,0.125,0,0.25,0.625,1,0.5,0.625,0,0.125,0,0.125,0.125,0.375,0.5,0.375,0,0.625,0.25,0.875,0.125,0.25 )
pB2 <- c(0,0.25,0.375,0.25,0.25,0.625,0.375,0.25,0,0.5,0.125,0.5,0.125,0.75,0.5,0.375,0.5,0.125,0,0.625,0.25,0.625,0,0.25,0.375 )

out <- cbind(A0,A1,A2,B0,B1,B2)
pro <- cbind(pA0,pA1,pA2,pB0,pB1,pB2)

max <- rep(-1,25)
min <- rep(9999,25)

max <- apply( cbind(pro,out) ,1, function(x) { 

	y <- cbind(x[1:6],x[7:12])
	z <- ifelse( y[,1] >0 , y[,2] , -99999 )
	max(z)
	
})

min <- apply( cbind(pro,out) ,1, function(x) { 

	y <- cbind(x[1:6],x[7:12])
	z <- ifelse( y[,1] >0 , y[,2] , 999999 )
	min(z)
	
})


# Define some population parameters, mean and standard devitions of r
# and fechner values
rm <- 0.65
rs <- 0.1

um <- 0.25
us <- 0.05

# Fechner will use a gamma distribution, so need to back out shape and
# scale parameters
k <- (um^2)/(us^2)
t <- (us^2)/um

# grab one value from this distribution
rval <-rnorm(1,mean = rm, sd = rs)
uval <-rgamma(1,shape = k, scale = t)

# R can have multiple datasets open at once. It calls these "data frames"
# I can think of a lot of instances where this would be hugely beneficial
d <- data.frame(pA0=pA0,pA1=pA1,pA2=pA2,pB0=pB0,pB1=pB1,pB2=pB2,
				A0=A0,A1=A1,A2=A2,B0=B0,B1=B1,B2=B2,
				max=max,min=min)

D <- d
D$ID <- 1
D$r  <-rnorm(1,mean = rm, sd = rs)
D$mu <-rgamma(1,shape = k, scale = t) 

for(i in 2:subjects){
	d$ID <- i
	d$r  <- rnorm(1,mean = rm, sd = rs)
	d$mu <- rgamma(1,shape = k, scale = t)

	D <- rbind(D,d)
}

# Practice using R as a functional language. This is a cheap example
# But you cannot do the equivaluent of this in Stata with "programs"
CRRA <- function(out,r){
	out^(1-r) / (1-r)
}

# Adding the Utility of A and B to the dataframe, this isn't actually
# necessary as it is in R. We could have just defined UA and UB as 
# vectors and they would have assumed the correct length. The benefit of
# not putting them in the dataframe would be that when the dataframe is
# saved, we don't keep these intermediary variables
D$ctx <- CRRA(max,D$r) - CRRA(min,D$r)

D$UA <- (pA0 * CRRA(A0,D$r) + pA1 * CRRA(A1,D$r) + pA2 * CRRA(A2,D$r)) / D$ctx / D$mu
D$UB <- (pB0 * CRRA(B0,D$r) + pB1 * CRRA(B1,D$r) + pB2 * CRRA(B2,D$r)) / D$ctx / D$mu

# Grab the probability of the choice of each row
D$pA <- exp(D$UA) / (exp(D$UA) + exp(D$UB))
D$pB <- exp(D$UB) / (exp(D$UA) + exp(D$UB))

# Random uniform number
rand <- runif(nrow(D))

# This is a great R function, ifelse collapses would would otherwise
# potentially be several lines of code into a very readable one line
# statement.
D$c <- ifelse(D$pA > rand, 0, 1)

save(D ,file="HO.Rda")

head(D,n=nrow(d))

pc<- ifelse(D$c==0,D$pA,D$pB)
mean(pc)

# Dataset is generated, call the file that does ML estimation.
#source("ML.r", echo=TRUE,print.eval=TRUE)
