# This is the main file to source the data and call the MSL functions

#libraries
library("parallel")

cores <- detectCores() 

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

mkhalton <- function(h=75,prime=3,HH=2,UH=1){

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

MSL <- function(par,h,h1,h2){

	par[2:4] <- exp(par[2:4])

	rm <- par[1]
	rs <- par[2]

	um <- par[3]
	us <- par[4]
	k <- (um^2)/(us^2)
	t <- (us^2)/um

	# Turn the halton sequences into distributions governed by par[]
	a <- matrix(qnorm(h1,mean=rm,sd=rs), ncol=h)
	b <- matrix(qgamma(h2,shape=k,scale=t), ncol=h) 

	dd <- data.frame( rbind(a,b) )

	sim <- lapply(dd,function(x){

		nn <- length(x) / 2
		r  <- 1 - x[1:nn]
		mu <- x[(nn+1):(2*nn)]

		ctx <- D$max^r/r - D$min^r/r

		UA <- (D$pA0 * D$A0^r/r + D$pA1 * D$A1^r/r) 
		UB <- (D$pB0 * D$B0^r/r + D$pB1 * D$B1^r/r) 

		# Calculte the probability stuff, need to rebase the utility because 
		# of precision and machine limit problems. exp(>709) evalutes as Inf.
		
		UB.1 <- (UB/ctx/mu) - (UA/ctx/mu)
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

    })

	# Each column of sim is a simulation, but we need the means of each row,
	# since each row corresponds to a unique choice. Thus the list needs to 
	# be made into a matrix to do sums. 

	like <- matrix(unlist(sim), ncol=h)
	sl <- sum(log(rowMeans(like)))

	return(-sl)

}

do.optim <- function(h,int){
    
	h1 <- H1[,1:h]
	h2 <- H2[,1:h]

	con <- list(trace=1,maxit=100)

	m <- optim(par=int,fn=MSL,h=h,h1=h1,h2=h2,method="BFGS", control=con,hessian=TRUE )

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
	pval<-2*(1-pt(abs(t),nrow(D)-length(int)))

	# Adjust for the logged parameters
	tt <-  m$par
	tt[2:4] <- exp(tt[2:4])

	ts <- se
	ts[2:4] <- exp(se[2:4])

	tu <- up
	tu[2:4] <- exp(tu[2:4])

	tl <- low
	tl[2:4] <- exp(tl[2:4])

	start <-int
	start[2:4] <- exp(int[2:4])

	# Save everything in a convienient place
	mm <- data.frame(real=real,init=start,est=m$par,par=tt,se=se,lower=tl,upper=tu,pvalue=pval,llike=m$value, H=h, HH=HH, UH=UH)
	# Print these things out
	print(mm)
	print("")
	return(mm)

}

# optim minimizes functions, so need to return the negitive
# of the log-likelihood in order to maximize the ll

real <- c(rm=mean(D$r),rs=sd(D$r),um=mean(D$mu),us=sd(D$mu))

# Which type of halton draws? 0 means every observation gets the same draw, 1
# means that each subject gets their own draw, 2 means every observaltion gets
# different draw.

allruns <- list()

for(HH in 0:2){

	UH <- 1
    
	H1 <- mkhalton(1000,prime=3,HH=HH,UH=UH)
	H2 <- mkhalton(1000,prime=7,HH=HH,UH=UH)

	models <- list()

	range <-list()
	jj <- 0
	for(i in seq(from=50,to=350,by=5)){
		jj <- jj + 1
		range[[jj]] <- i
	}

	dist <- real
	binit <- dist
	binit[2:4] <- log(dist[2:4])

	bcs <- mclapply(range,do.optim,int=binit,mc.cores=cores )

	models[[1]] <- bcs

	dist <- c(rm=.5,rs=.5,um=.5,us=.5)
	ainit <- dist
	ainit[2:4] <- log(dist[2:4])

	a50 <- mclapply(range,do.optim,int=ainit,mc.cores=cores )

	models[[2]] <- a50

	dist <- c(rm=runif(1),rs=runif(1),um=runif(1,min=.2),us=runif(1))
	uinit <- dist
	uinit[2:4] <- log(dist[2:4])

	runi <- mclapply(range,do.optim,int=uinit,mc.cores=cores )

	models[[3]] <- runi

	# And now for the same ones with a uniform genorator

	UH <- 0

	H1 <- mkhalton(1000,prime=3,HH=HH,UH=UH)
	H2 <- mkhalton(1000,prime=7,HH=HH,UH=UH)

	bcs <- mclapply(range,do.optim,int=binit,mc.cores=cores )

	models[[4]] <- bcs

	a50 <- mclapply(range,do.optim,int=ainit,mc.cores=cores )

	models[[5]] <- a50

	runi <- mclapply(range,do.optim,int=uinit,mc.cores=cores )

	models[[6]] <- runi

	# Done with this H

	allruns[[HH+1]] <- models

}

save(allruns,file="allruns.Rda")

