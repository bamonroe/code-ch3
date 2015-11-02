# Draw up graphs in R to demonstrate how the RP model operates with the example
# given in Ch.2, The Stochastic Money Pump

# Load in some libraries:
# rootSolve finds all 0s of a nonlinear equation
library(rootSolve)
# Let's try and do this with a modern graphics package
library(ggplot2)

# The instrument, framed in perspective of the subject, thus the buy amount is
# the amount that the subject will pay to buy the lottery.

pA0 <- .5
pA1 <- .5

pB0 <- 1

A0 <- 10
A1 <- 100

Buy  <- 55.5
Sell <- 54.5

# Utility functions
crra <- function(x,r){
    ifelse(r==1, log(x),  x^(1-r) / (1-r) )
 } 
UA <- function(r){
	(pA0 * crra(A0,r) + pA1 * crra(A1,r)) 

}

UBuy <- function(r){
	crra(Buy,r)
}

# Certainty Equivalents

CE <- function(u,r){
	ifelse(r==1 , exp(u) ,(u * (1-r) ) ^ (1 / (1-r)) )
}

CEA <- function(r){
	CE(UA(r),r)
}
CEBuy <- function(r){
	CE(UBuy(r),r)
}

CESell <- function(r){
	CE(USell(r),r)
}


# Some sensible boundaries for the plot
ybot <- -55
ytop <- 85

xbot <- -.1
xtop <- 1.65


# make a dataset to plot the data - didn't want to do it this way, but I can't
# find an alternative way to get a polygon underneath the the curves...

domain <- data.frame(r=seq(from=xbot, to=xtop, by=.0001))
attach(domain)
domain$UA   <- UA(r)
domain$UBuy <- UBuy(r)
detach(domain)

# A function which will return 0 if the utilities of two options, given two,
# presumably different, rvalues, are equal
solve <- function(F1,F2,r1,r2=r1){
	F1(r1) - F2(r2)
}

# A vector values to plot. These are the rvalues associated with the draw
# for the utility of the ticket. 

rvals <- c(.5)

# what are the rvals of area to shade?

shade0 <- subset(domain , domain$UBuy < UA(rvals[1]) & domain$r < 1  )

shade1 <- subset(domain , domain$UBuy < UA(rvals[1]) & domain$r > 1  )
shade1 <- rbind( 
				shade1 ,  
				c(shade1$r[nrow(shade1)],UA(rvals[1]),UBuy(rvals[1])) , 
				c(shade1$r[1],UA(rvals[1]),UBuy(rvals[1])) , 
				c(shade1$r[1],ybot,ybot)
				)

if(nrow(shade0) > 0){
#	shade0 <- rbind(
#				c(xbot,shade0$UA[nrow(shade0)],shade0$UB[nrow(shade0)]),
#				shade0
#			)
	shade <- rbind(shade0,shade1)
} else{
	shade <- shade1
}

# Make the first element of rvals the r value which indicates indifference
x0 <- rvals

y0 <- rep(ybot, length(x0))

x1 <- c( x0, x0 )
y1 <- c(y0 ,UA(x0))

# Make every observation in the dataset a line, for every rval, there are 2
# verticle lines and 1 horizontal line

l0 <-	sapply( x0,function(x) {  
			uniroot.all(solve,c(-1,1),F1=UA,F2=UA,r1=x)[2] 
		})

x1 <- c(x1 , l0)
y1 <- c(y1 , y0)

x2 <- c(x0 , l0 , l0 )
y2 <- c(UA(x0))

lines <- data.frame(x1=x1 , x2=x2 , y1=y1 , y2=y2)

# Some calculations for the area which is 0 probability
minBuy <- optim(par=.7,UBuy,method="Brent",lower=0,upper=1)
root0 <- uniroot.all(solve,c(-10,1),F1=UA,F2=UBuy,r2=minBuy$par)
nl0  <- c(root0[1],root0[1],ybot,UA(root0[1]))
nl1  <- c(root0[2],root0[2],ybot,UA(root0[2]))
nl2  <- c(root0[1],root0[2],UA(root0[1]),UA(root0[2]))

lines <- rbind(lines,nl0,nl1,nl2)

#Make the sekelton 
p <- ggplot(data = domain, mapping = aes(x=x)) 

# Add the polygons
p <- p + geom_polygon(data=shade0, aes(x=r, y=UBuy ),fill="dark green", alpha=.51)
p <- p + geom_polygon(data=shade1, aes(x=r, y=UA ),fill="dark green", alpha=.51)

# Plot the lines
p <- p + geom_segment(data = lines ,aes(x = x1, y = y1, xend = x2, yend = y2 ) , linetype="longdash")  

#plot the functions
p <- p + geom_line( aes(x = r, y = UA, color = "blue"))
p <- p + geom_line( aes(x = r, y = UBuy, color = "red"))

# Plot r == 1
points <- subset(domain , domain$r == 1  )
p <- p + geom_point( data=points, aes(x = r, y = UA, color = "blue"))
p <- p + geom_point( data=points, aes(x = r, y = UBuy, color = "red"))

# Make better looking

p <- p + scale_y_continuous(
    				   limits = c(ybot, ytop),
    				   breaks=c(seq(from=ybot, to=ytop, by=15))
    				   ) +
    scale_x_continuous(
    				   limits = c(xbot, xtop),
    				   breaks=c(seq(from=xbot, to=xtop ,by=.1),1)
    				   ) +
    scale_color_manual(
    				   name = "Functions",
                       values = c("blue", "red"), # Color specification
                       labels = c("Utility of Ticket", "Utility of Buy Price","Utility of Sell Price")
                       ) +
	theme(
		axis.text.x=element_text(angle=50, size=12, vjust=0.5), 
		axis.title.x = element_text(color="forestgreen", vjust=-0.35),
		axis.title.y = element_text(color="forestgreen" , vjust=0.35)
		) +
	labs(x="CRRA Value", y="Utility", title="Buying the Ticket") 

p
ggsave(file="graph.png")
ggsave(file="graph.pdf")
