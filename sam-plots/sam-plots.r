# Clear All
rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(parallel)

cores <- detectCores()

# Grab our data - A dataframe called "D"
load("../data/sam-dat/2.5Mil.Rda")

# Transform it into a table dataframe for use with dplyr
D <- tbl_df(D)

# Need to make at least 4 groups
#	- Consistent
#	- FOSD Only
#	- Light MSB
#	- Light MSB + FOSD

D <- arrange(D,Pattern)

# Recreate the pattern variable in D, this actually is backwards from the way it is in D,
# but it is still the most efficient way that I'm aware of.
patmat <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
pats <- data.frame(t(patmat))
pat2 <- patmat

# Just flip the thing around
for(i in 10:1){
	pat2[,i] <- patmat[,11-i]
}
patmat <- pat2
rm(pat2)


getCategory <- function(pat){

	cat <- c()

	# FOSD is easy, 0 in row 10 is FOSD
	FOSD <- ifelse(pat[10]==0,T,F)

	# Consistent means 1 switch point with 0 in row 1 or 0 switch points
	switches <- 0

	for(i in 2:9){
		if(pat[i] != pat[i-1]) switches <- switches + 1
	}

	consistent <- (pat[1] == 0 & switches ==1 ) | (pat[1] ==0 & switches ==0 ) | (pat[1] ==1 & switches ==0 )

	LightMSB <- (pat[1]==0 & switches == 3 ) | (pat[1]==0 & switches == 2 )

	if( consistent & !FOSD ){
		cat <- c(0)
	}
	else if( consistent & FOSD ){
		cat <- c(1)
	}
	else if( LightMSB & !FOSD ){
		cat <- c(2)
	}
	else if( LightMSB & FOSD ){
		cat <- c(3)
	}
	else{
		cat <- c(4)
	}

	cat

}

D$Type <- apply(patmat,1,getCategory)
E <- filter(D,Type < 4)

# Use arrange so that Consistent is plotted last
E <- arrange(E,desc(Type))

E$Type <- factor(E$Type, labels=c("Consistent","FOSD Only","Light MSB","Light MSB + FOSD"))

E$ll <- log(E$PC)

# There are three points to make, A, B and C. A is the most likely consistent point
# B is the most likely FOSD only point, C is the consistent point that is closest in
# liklihood to B that has greater likelihood

cons <- E %>%
	select(Type,WP,ll) %>%
	filter(Type=="Consistent") %>%
	arrange(desc(ll))

FOSD <- E %>%
	select(Type,WP,ll) %>%
	filter(Type=="FOSD Only") %>%
	arrange(desc(ll))

# A and B are straight forward
A.x <- cons$ll[1]
A.y <- cons$WP[1]

B.x <- FOSD$ll[1]
B.y <- FOSD$WP[1]

# Need to find the 

diff <- cons %>%
	mutate(diff =(ll - B.x)) %>%
	filter(diff>0) %>%
	arrange(diff)

C.x <- diff$ll[1]
C.y <- diff$WP[1]



stop("here")



# Now to configure the graph aesthetically

# Configure the points
## Shape of the points
shape <- 20
## Size of the points
point.size <- 5
## alpha
alph <- .9

# Configure x-axis title
## title size
x.title.size <- 14
## title vertical justification
x.title.vjust <- 1

# Configure x-axis text
## text size
x.text.size <- 11
## angle of the text
x.text.angle <- 45


# Configure y-axis title
## title size
y.title.size <- 14
## title horizontal justification
y.title.hjust <- .5
## title vertical justification
y.title.vjust <- 1

# Configure y-axis text
## text size
y.text.size <- 11
## angle of the text
y.text.angle <- 45







p <- ggplot(data=E, aes_string(x="ll"))
p <- p + geom_point(shape=shape,size=point.size, alpha=alph, aes_string(y="WP",color="Type") )
p <- p + labs(x="Log of Simulated Likelihood",y="Expected Ratio of Obtained to Optimal Welfare")
p <- p + theme(axis.title.x=element_text(size=x.title.size),
			   axis.text.x=element_text(size=x.text.size,angle=x.text.angle),
			   axis.title.y=element_text(size=y.title.size,hjust=y.title.hjust,vjust=y.title.vjust),
			   axis.text.y=element_text(size=y.text.size,angle=y.text.angle)
			)

p

stop("here")

# Colors
cy  <-"#32b0e6" # Cyan
yo  <-"#e6c832" # Yellow Orange
ma  <-"#e63271" # Magenta

lcy <- "#41c6ff"
dyo <- "#ff7d41"
ye  <- "#ffdc41"
pur <- "#b941ff"

neog <- "#7aff32"
mpin <- "#ff327a"
lora <- "#ff7a32"
skyb <- "#32beff"

# 3 tone colors
color.3 <- c(cy,yo,ma)
color.3.1 <- c("red","blue","yellow")
# 4 tone colors
color.4.0 <- c(lcy,dyo,pur,ye)
color.4.1 <- c(ye,lcy,dyo,pur)

color.4.2 <- c(cy,ye,dyo,pur)
color.4.3 <- c(mpin,neog,lora,skyb)

color.4.4 <- c(dyo,cy,yo,pur)

# Colors to use
colors <- color.4.4

# Shape of the points
shape <- 20

# Faceting, 
USE.w <- melt(get(USE), measure.vars=c("M.WP","V.WP","M.WC","V.WC") )
USE.e <- melt(get(USE), measure.vars=c("M.EE","V.EE","M.WE0","V.WE0") )

USE.w$variable  <-  ifelse(USE.w$variable == "M.WP" , "Mean Expected Welfare Proportion",
					ifelse(USE.w$variable == "V.WP" , "Variance of Expected Welfare Proportion",
					ifelse(USE.w$variable == "M.WC" , "Mean Expected Welfare Surplus",
													  "Variance of Expected Welfare Surplus")
					))

USE.e$variable  <-  ifelse(USE.e$variable == "M.EE" , "Mean Expected Number of Errors",
					ifelse(USE.e$variable == "V.EE" , "Variance of Expected Number of Errors",
					ifelse(USE.e$variable == "M.WE0", "Mean Expected Proportion of No Error Choices",
													  "Variance of Expected Proportion of No Error Choices")
					))


USE.w$variable <- factor(USE.w$variable, levels=c(	"Mean Expected Welfare Proportion",                              											  
													"Variance of Expected Welfare Proportion",
													"Mean Expected Welfare Surplus",
													"Variance of Expected Welfare Surplus"))

USE.e$variable <- factor(USE.e$variable, levels=c(	"Mean Expected Number of Errors",
													"Variance of Expected Number of Errors",
													"Mean Expected Proportion of No Error Choices",
													"Variance of Expected Proportion of No Error Choices"))
													

# What are the parameter names
names   <- c("rm","rs","um","us")
s.names <- c("rs","rm","us","um")

names   <- rep(names,2)
s.names <- rep(s.names,2)

p.t <- c(rep("USE.e",4),rep("USE.w",4))

to.plot <- data.frame(rbind(names,s.names,p.t))

getPlotted <- function(plot){
    
    x.par <- as.character(plot[1])
    s.par <- as.character(plot[2])
    data  <- as.character(plot[3])

	# Reference point for scale
	mid.ref <- (max(get(data)[[s.par]]) - ((max(get(data)[[s.par]])-(min(get(data)[[s.par]])))/2))
	limits <- c(min(get(data)[[s.par]]), max(get(data)[[s.par]]))

	q.dist <- (limits[2] - mid.ref) / 2

	bbr <- mid.ref - q.dist
	tbr <- mid.ref + q.dist

	breaks <- c(bbr,mid.ref,tbr,limits[2])
	breaks <- floor(breaks * 100) / 100
	ll <- ceiling(limits[1]*100)/100
	breaks <- c(ll,breaks)


	# Values for the various plots
	title <- "This is a title\n\n\n"

	leg.title.size  <- 22
	leg.title.vjust <- -100

	leg.text.size   <- 14
	leg.text.angle  <- 45
	leg.position    <- "right"
	leg.direction	<- "horizontal"
	leg.direction	<- "vertical"
	leg.key.height  <- 1
	leg.key.width   <- .25

	x.title <- ifelse(x.par=="rm", "Mean of CRRA",
			   ifelse(x.par=="rs", "Standard Deviation of CRRA",
			   ifelse(x.par=="um", "Mean of Lambda","Standard Deviation of Lamda")))

	x.title.size <- 24


	leg.title <- ifelse(s.par=="rm", "Mean of CRRA",
                 ifelse(s.par=="rs", "Standard Deviation of CRRA",
                 ifelse(s.par=="um", "Mean of Lambda","Standard Deviation of Lamda")))

	fac.size <- 24

	# gather the plot layers
	p <- ggplot(data=get(data), aes_string(x=x.par))
	p <- p + scale_color_gradientn(name=leg.title,space="Lab",colours=colors,limits=limits,breaks=breaks) + 
	  	      scale_fill_gradientn(name=leg.title,space="Lab",colours=colors,limits=limits,breaks=breaks)
	p <- p + geom_point(shape=shape, alpha=alph, aes_string(y="value",color=s.par,fill=s.par) )
	p <- p + facet_wrap( facets=~variable, ncol=2, scale="free_y")
	p <- p + labs(x=x.title,y=NULL)
	p <- p + theme(legend.title=element_text(size=leg.title.size,vjust=.9),
					legend.text=element_text(size=leg.text.size,angle=leg.text.angle),
					legend.position=leg.position,
					legend.direction=leg.direction,
					legend.key.height=unit(leg.key.height,"in"),
					legend.key.width=unit(leg.key.width,"in"),
					axis.title.x=element_text(size=x.title.size),
					strip.text=element_text(size=fac.size)
					)

	# Return the plot to the list
	dat <- ifelse(data=="USE.w","Wel","Err")

	name <- paste(x.par,"-",dat,sep="")

	ret <- c(name=list(x.par,dat,p))

	# Diminetions of plots
	h <- 8.50 - (.79*2) - .5	# Page is 8.5 x 11 inches, margins are .79 inches
	w <- 11.0 - (.79*2)

	save.scale <- 3

	fname <- paste("../data/agg-plots/",dat,"-",x.par,".jpg",sep="")

	ggsave(filename=fname, plot=p, width=w, height=h, units="in",scale=save.scale )

	return(ret)

}

plots <- mclapply(to.plot,getPlotted,mc.cores=cores)
#plots <- lapply(to.plot,getPlotted)


