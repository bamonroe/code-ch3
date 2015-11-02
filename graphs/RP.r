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
UC <- function(r){
	(pA0 * crra(20,r) + pA1 * crra(80,r)) 
}

UBuy <- function(r){
	crra(Buy,r)
}

USell <- function(r){
	crra(Sell,r) 
}

# Certainty Equivalents

CE <- function(u,r){
	(u * (1-r) ) ^ (1 / (1-r))
}

CEA <- function(r){
	CE(UA(r),r)
}
CEC <- function(r){
	CE(UC(r),r)
}
CEBuy <- function(r){
	CE(UBuy(r),r)
}

CESell <- function(r){
	CE(USell(r),r)
}

# Some sensible boundaries for the plot
ybot <- 0
ytop <- 85

# A function which will return 0 if the utilities of two options, given two,
# presumably different, rvalues
solve <- function(F1,F2,r1,r2=r1){
	F1(r1) - F2(r2)
}

#A function to be used to get the roots of values
groots <- function(r){
	uniroot.all(solve,c(-10,1),F1=UBuy,F2=UA,r2=r)
}

# A vector values to plot. These are the rvalues associated with the draw
# for the utility of the ticket. They should all be less than the value which
# makes the minimum for UA ~ .75

# Make the first element of rvals the r value which indicates indifference
x0 <- c( uniroot.all(solve,c(-10,1),F1=UBuy,F2=UA), seq(from=0.05, to = .5, by = .2 ))
#x0 <- c( uniroot.all(solve,c(-10,1),F1=UBuy,F2=UA))
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

# Start filling in the polygons
# When using geom_polygon, you will typically need two data frames:
# one contains the coordinates of each polygon (positions),  and the
# other the values associated with each polygon (values).  An id
# variable links the two together

# lets see if we can do this a bit more simply and more generalized
# We have only the horizontal lines from before, x1 will be the x cordinate for
# the left side of the polygon and x2 will be for the right. 

n <- nrow(lines[which(lines$x1 != lines$x2),])

hl <- data.frame(
	x = c(
		lines$x1[which(lines$x1 != lines$x2)],
		lines$x1[which(lines$x1 != lines$x2)],
		lines$x2[which(lines$x1 != lines$x2)],
		lines$x2[which(lines$x1 != lines$x2)]
		),
	y = c(rep(ybot,n),rep(ytop,n),rep(ytop,n),rep(ybot,n)),
	id =c(1:n)
)

#Make the sekelton 
p <- ggplot(data = data.frame(x = 0), mapping = aes(x=x)) 

# Add the polygons
p <- p + geom_polygon(data=hl, aes(x=x, y=y,group=id ), fill="dark green",alpha=.15)

# Plot the lines
p <- p + geom_segment(data = lines ,aes(x = x1, y = y1, xend = x2, yend = y2 ) , linetype="longdash")  

# Put the functions in
p <- p +
	layer(stat = "function",
          fun = UA,
          mapping = aes(color = "UA")
          ) +
    layer(stat = "function",
          fun = UBuy,
          mapping = aes(color = "UBuy")
          ) +
    scale_x_continuous(
    				   limits = c(-.25, 1),
    				   breaks=c(seq(from=-.2, to=1, by=.15))
    				   ) +
    scale_y_continuous(
    				   limits = c(ybot, ytop),
    				   breaks=c(seq(from=ybot, to=ytop, by=15))
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

#Plot the plot
p
