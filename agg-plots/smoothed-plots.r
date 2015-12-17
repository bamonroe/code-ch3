# Clear All
rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(grid)
library(reshape2)
library(parallel)

cores <- detectCores()

## Haven't use, but will likely use
library(np)

# Grab our data
load("../data/agg-dat/Agg1Mil-S10k.Rda")
MM <- tbl_df(MM)

# Don't always want to use it all, there is tons of data
sam.prop <- .21

# Bounds for the means of data
lbound <- -1.7134
ubound <- 1.3684

lbound <- -1.9
ubound <- 1.55

# Certain data have rm inside the relevent range of the HL-MPL, others don't.
# Make the distiction
IN  <- MM %>%
	sample_frac(sam.prop) %>%
	filter(rm > lbound, rm < ubound) %>%
	mutate(id = 1:n())

OUT <- MM %>%
	sample_frac(sam.prop) %>%
	filter(rm < lbound | rm > ubound) %>%
	mutate(id = 1:n())

# Functionally determine the alpha channel for the points
alph <- (.75 / (exp(sam.prop))) 
#alph <- .35
	
# Text naming the data.frame to use
USE <- "IN"

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
    data  <- get(as.character(plot[3]))

	# Reference point for scale
	mid.ref <- (max(data[[s.par]]) - ((max(data[[s.par]])-(min(data[[s.par]])))/2))
	limits <- c(min(data[[s.par]]), max(data[[s.par]]))

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

	# Split the dataset up in a few ways to get multiple smoothed lines
	s.num <- 5

	data <- data %>% 
			mutate_(.dots=setNames(paste0("ntile(",s.par,",",s.num,")"),"bin"))

	splits <- dlply(data,"bin")

	p <- qplot()
	for( i in 1:s.num){
		#p <- p + geom_smooth(data=splits[[i]],stat="smooth", method="gam", aes_string(x=x.par,y="value"), span=0.01 , formula=y~s(x^3) ) 
		p <- p + stat_smooth(data=splits[[i]],geom="smooth", method="gam", aes_string(x=x.par,y="value"), span=0.01 , formula=y~s(x^3) ,color=(i)) 
	}
	p <- p + facet_wrap( facets=~variable, ncol=2, scale="free_y")

#	# gather the plot layers
#	p <- ggplot(data=data, aes_string(x=x.par))
#	p <- p + scale_color_gradientn(name=leg.title,space="Lab",colours=colors,limits=limits,breaks=breaks) + 
#	  	      scale_fill_gradientn(name=leg.title,space="Lab",colours=colors,limits=limits,breaks=breaks)
#	#p <- p + geom_point(shape=shape, alpha=alph, aes_string(y="value",color=s.par,fill=s.par) )
#
#	for( i in 1:s.num){
#		p <- p + geom_smooth(data=splits[[i]],stat="smooth", method="gam", aes_string(x=x.par,y="value"), span=0.9 , formula=y~s(x^2) ) 
#	}
#
#	#p <- p + geom_smooth(stat="smooth", span=.1, aes_string(y="value",color="bin") )
#
#	p <- p + facet_wrap( facets=~variable, ncol=2, scale="free_y")
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
	dat <- ifelse(as.character(plot[3])=="USE.w","Wel","Err")

	name <- paste(x.par,"-",dat,sep="")

	ret <- c(name=list(x.par,dat,p))

	# Diminetions of plots
	h <- 8.50 - (.79*2) - .5	# Page is 8.5 x 11 inches, margins are .79 inches
	w <- 11.0 - (.79*2)

	save.scale <- 3

	fname <- paste("../data/agg-plots/Smooth-",dat,"-",x.par,".jpg",sep="")

#	ggsave(filename=fname, plot=p, width=w, height=h, units="in",scale=save.scale )

	return(ret)

}

#plots <- mclapply(to.plot,getPlotted,mc.cores=cores)
plots <- lapply(to.plot,getPlotted)

plots[[1]]

