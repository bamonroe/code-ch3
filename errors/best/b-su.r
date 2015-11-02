# Recreate Error Calculation
# Clear all
rm(list = ls())

GetRes <- function(sim){

	# The passed r value is the real "r", but all the calculations use 1-r
	# whenever r is needed, so subtract it once initially and save the 
	# computation times of all the subsequent calls, these tiny fractions
	# of a second add up with a very large snum.
	r1  <- 1 - sim[1]
	mu <- sim[2]

	# Calculate the basline utility
	UA <- (pA0 * A0^r1/r1 + pA1 * A1^r1/r1)
	UB <- (pB0 * B0^r1/r1 + pB1 * B1^r1/r1)

	# Figure out which of these are errors
	c.AB <- UA > UB & cset == 1
	c.BA <- UB > UA & cset == 0
	e1 <- c.AB | c.BA

	c.err <- ifelse(e1,1,0)

	# Make some vectors of certainty equivalents
	CEA <- (UA * r1)^(1/r1)
	CEB <- (UB * r1)^(1/r1)

	# Get the CE of the chosen option
	W0 <- ifelse(cset==0,CEA,CEB)

	# Get th CE of the unchosen option
	W1 <- ifelse(cset==1,CEA,CEB)

	# Get the CE of the Error only
	WE <- ifelse(c.err==1,W0-W1,0)

	# Get the maximum CE
	WMax <- ifelse(CEA>CEB,CEA,CEB)

	# Calculte the probability stuff, need to rebase the utility because 
	# of precision and machine limit problems. exp(>709) evalutes as Inf.
	
	UB.1 <- UB/mu - UA/mu
	UA.1 <- 0

	# Have things gone haywire because of the computer's inability to handle
	# numbers bigger than ~3e310 ?
	c.N <- is.nan(UB.1) | is.infinite(UB.1)

	pA <- ifelse( c.N ,		# Are we dealing with an insane number?
		# yes:
		ifelse( UB > UA , 0 , 1 ) ,
		# no, but are we making an insane number via exp?
		ifelse( UB.1 > 709 , 0 , exp(UA.1) / (exp(UA.1) + exp(UB.1)) )
	)

	pB <- 1 - pA

	# Store the probabilities of the chosen options
	p0 <- ifelse(cset==0,pA,pB)

    # Gather the results 
    res <- c(0,0,0,0,0)

	# How many errors
    res[1] <- sum(c.err)

    # Probability of this pattern
    res[2] <- prod(p0)/snum

    # Welfare metric mirroring opportunity cost of whole choice set
    res[3] <- sum(W0 - W1)/snum

    # Welfare metric for proportion of welfare obtained / total availavle
	res[4] <- sum(W0)/sum(WMax)/snum

	# Welfare metric opportunity cost of error only
	res[5] <- sum(WE)/snum

	return(res)

	# 1- Err# , 2-Prob, 3-WOC, 4-Wprop
    
}

mkpattern <- function(){

	pattern <- matrix(nrow = 1024, ncol = 11) 
	colnames(pattern) <- list("count","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10")

	count <- 0
    
	for(i1 in 0:1){
		for(i2 in 0:1){
			for(i3 in 0:1){
				for(i4 in 0:1){
					for(i5 in 0:1){
						for(i6 in 0:1){
							for(i7 in 0:1){
								for(i8 in 0:1){
									for(i9 in 0:1){
										for(i10 in 0:1){

											count <- count + 1
											pattern[count][1] <- count

											pattern[count,2] <- i1
											pattern[count,3] <- i2
											pattern[count,4] <- i3
											pattern[count,5] <- i4
											pattern[count,6] <- i5
											pattern[count,7] <- i6
											pattern[count,8] <- i7
											pattern[count,9] <- i8
											pattern[count,10] <- i9
											pattern[count,11] <- i10

										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return(pattern)

}

pattern <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))

# Set up our instrument
A0 <- rep(1.60,10)
A1 <- rep(2.00,10)
B0 <- rep(0.10,10)
B1 <- rep(3.85,10)

pA0 <- seq(.9,0,by=-.1)
pB0 <- seq(.9,0,by=-.1)
pA1 <- seq(.1,1,by=.1)
pB1 <- seq(.1,1,by=.1)

Max <- rep(3.85,10)
Min <- c(rep(0.10,9),2)

# What is the number of simulations to run
snum <- 100
# How many choices per subject
cnum <- length(A0)

# Keep track of which simulation we're on
dnum <- 0

# How many iterations to run?
iternum <- 1

# There are 15 different versions I want to run. These hold groups of
# parameters constant. This will hopefully add some statistical power

vnum <- 1

#for(iter in 1:iternum){
while(TRUE){
    
	start <- proc.time()
    
	dnum <- dnum + 1
	
	# Set up a DataFram to hold all the error data
	D <- data.frame(E.0=rep(0,1024),E.1=rep(0,1024),E.2=rep(0,1024),E.3=rep(0,1024),E.4=rep(0,1024),
		E.5=rep(0,1024),E.6=rep(0,1024),E.7=rep(0,1024),E.8=rep(0,1024),E.9=rep(0,1024),E.10=rep(0,1024)
		)

	D$EE <- 0
	D$PC <- 0
	D$WP <- 0
	D$WC <- 0
	D$WE <- 0

	# For every simulation, I want the number of errors. I also want this per
	# pattern, so every subject will get a 

	# Define some population parameters, mean and standard devitions of r
	# and fechner values

	# I'm thinking it best to use the same distribution of r with different distributions of fechner
	# to increase statistical power. This can be a pain if I take this approach, so each r distribution
	# will only be repeated once. Will do the same with fechner, but with a different file.

	print(paste("We're on dataset #",dnum," and vnum #",vnum,sep=""))

	rm.min <- -1.75
	rm.max <- 1.75
	um.min <- .05
	um.max <- 3.5
	rs.min <- .01
	rs.max <- 3.5
	us.min <- .01
	us.max <- 4

	if(dnum %% 2 == 1){
		# Everytime we have an odd number, make every parameter a fresh one
		rm <- runif(1,min=rm.min,max=rm.max)
		rs <- runif(1,min=rs.min,max=rs.max)
		um <- runif(1,min=um.min,max=um.max)
		us <- runif(1,min=us.min,max=us.max)
	}
	else{
		# Everytime we have an even number, only change certain parameters
		if(vnum == 1){
#			rm <- runif(1,min=rm.min,max=rm.max)
			rs <- runif(1,min=rs.min,max=rs.max)
			um <- runif(1,min=um.min,max=um.max)
			us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 2){
	#		rm <- runif(1,min=rm.min,max=rm.max)
	#		rs <- runif(1,min=rs.min,max=rs.max)
			um <- runif(1,min=um.min,max=um.max)
			us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 3){
	#		rm <- runif(1,min=rm.min,max=rm.max)
	#		rs <- runif(1,min=rs.min,max=rs.max)
	#		um <- runif(1,min=um.min,max=um.max)
			us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 4){
	#		rm <- runif(1,min=rm.min,max=rm.max)
	#		rs <- runif(1,min=rs.min,max=rs.max)
			um <- runif(1,min=um.min,max=um.max)
	#		us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 5){
	#		rm <- runif(1,min=rm.min,max=rm.max)
			rs <- runif(1,min=rs.min,max=rs.max)
	#		um <- runif(1,min=um.min,max=um.max)
			us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 6){
	#		rm <- runif(1,min=rm.min,max=rm.max)
			rs <- runif(1,min=rs.min,max=rs.max)
	#		um <- runif(1,min=um.min,max=um.max)
	#		us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 7){
	#		rm <- runif(1,min=rm.min,max=rm.max)
			rs <- runif(1,min=rs.min,max=rs.max)
			um <- runif(1,min=um.min,max=um.max)
	#		us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 8){
			rm <- runif(1,min=rm.min,max=rm.max)
	#		rs <- runif(1,min=rs.min,max=rs.max)
			um <- runif(1,min=um.min,max=um.max)
			us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 9){
			rm <- runif(1,min=rm.min,max=rm.max)
	#		rs <- runif(1,min=rs.min,max=rs.max)
	#		um <- runif(1,min=um.min,max=um.max)
			us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 10){
			rm <- runif(1,min=rm.min,max=rm.max)
	#		rs <- runif(1,min=rs.min,max=rs.max)
			um <- runif(1,min=um.min,max=um.max)
	#		us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 11){
			rm <- runif(1,min=rm.min,max=rm.max)
	#		rs <- runif(1,min=rs.min,max=rs.max)
	#		um <- runif(1,min=um.min,max=um.max)
	#		us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 12){
			rm <- runif(1,min=rm.min,max=rm.max)
			rs <- runif(1,min=rs.min,max=rs.max)
			um <- runif(1,min=um.min,max=um.max)
	#		us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 13){
			rm <- runif(1,min=rm.min,max=rm.max)
			rs <- runif(1,min=rs.min,max=rs.max)
	#		um <- runif(1,min=um.min,max=um.max)
	#		us <- runif(1,min=us.min,max=us.max)
		}

		else if(vnum == 14){
			rm <- runif(1,min=rm.min,max=rm.max)
			rs <- runif(1,min=rs.min,max=rs.max)
			um <- runif(1,min=um.min,max=um.max)
	#		us <- runif(1,min=us.min,max=us.max)
		}

		else if (vnum == 15){
			rm <- runif(1,min=rm.min,max=rm.max)
			rs <- runif(1,min=rs.min,max=rs.max)
			um <- runif(1,min=um.min,max=um.max)
			us <- runif(1,min=us.min,max=us.max)
			vnum <- 0
		}
		vnum <- vnum + 1
	}
	
	D$rmean  <- rm
	D$mumean <- um
	D$rstd   <- rs
	D$mustd  <- us
	D$ID     <- pattern[,1]
	D$pattern <-  apply(pattern,1,toString)

	# Fechner will use a gamma distribution, so need to back out shape and
	# scale parameters
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	sim <- matrix(  c(rnorm(snum,mean = rm, sd = rs), rgamma(snum,shape = k, scale = t) ) , nrow=snum  )

	# Give the very wide bounds we've mad for ourselves, occasionally we get 
	# a fechner value equal to 0. This creates an "oh shi-" problem later on.
	# Keep replacing 0s with another random draw from the same distribution until
	# they go away.
	while(any(sim[,2]==0)){
		sim[,2] <- ifelse(sim[,2]==0 , rgamma(1,shape = k, scale = t),sim[,2] )
	}

	for(i in 1:nrow(pattern)){

		print(i)

		cset <- matrix(pattern[i,],nrow=10)

		# 1- Err# , 2-Prob, 3-WOC, 4-Wprop

		result <- apply(sim,1,GetRes)

		D$EE[i] <- sum(result[1,])/snum
		D$PC[i] <- sum(result[2,])

		D$WC[i] <- sum(result[3,])
		D$WP[i] <- sum(result[4,])
		D$WE[i] <- sum(result[5,])

		D$E.0[i] <- sum(result[1,] == 0)
		D$E.1[i] <- sum(result[1,] == 1)
		D$E.2[i] <- sum(result[1,] == 2)
		D$E.3[i] <- sum(result[1,] == 3)
		D$E.4[i] <- sum(result[1,] == 4)
		D$E.5[i] <- sum(result[1,] == 5)
		D$E.6[i] <- sum(result[1,] == 6)
		D$E.7[i] <- sum(result[1,] == 7)
		D$E.8[i] <- sum(result[1,] == 8)
		D$E.9[i] <- sum(result[1,] == 9)
		D$E.10[i] <-sum(result[1,] == 10)

	}

	save(D ,file=paste("Dat",dnum,"-SU.Rda",sep=""))
	rm(D)

	end <- proc.time()
	time <- end - start
	print(time)

	conn <- file("../breaker.txt", "rt")
	bk <- substr(scan(conn, what = "", sep = "\n", quiet = TRUE), 1, 1)
	close(conn)
	if(	bk == "1")	break

}

