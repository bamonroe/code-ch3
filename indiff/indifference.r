# Get the indifference points of lotteries
#Using the Holt Laury(2002) instrument
A0 <- rep(1.60,10)
A1 <- rep(2.00,10)
B0 <- rep(0.10,10)
B1 <- rep(3.85,10)

pA0 <- seq(from=0.9, to=0, by= -.1)
pB0 <- seq(from=0.9, to=0, by= -.1)
pA1 <- seq(from=0.1, to=1, by= .1)
pB1 <- seq(from=0.1, to=1, by= .1)

inst <- matrix(c(A0,A1,B0,B1,pA0,pA1,pB0,pB1),nrow=10)

mindiff <- function(r,lots){
	a0  <- lots[1]
	a1  <- lots[2]
	b0  <- lots[3]
	b1  <- lots[4]
	pa0 <- lots[5]
	pa1 <- lots[6]
	pb0 <- lots[7]
	pb1 <- lots[8]

	UA <- (a0^(1-r))/(1-r)*pa0 + (a1^(1-r))/(1-r)*pa1
	UB <- (b0^(1-r))/(1-r)*pb0 + (b1^(1-r))/(1-r)*pb1

	diff <- (UA -UB)^2

	return(diff)
}

indiff <- apply(inst,1,function(x){
					
		print(x)
		bound <- optimize(f = mindiff,interval=c(-2,2),lots=x)
		bound$minimum

})

