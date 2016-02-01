# Clear All
rm(list=ls())

use.par <- F

if(use.par){
	# Add in Rhpc support functions
	source("../Rhpc/cl-tools.r")
}

add.lib <- function(){
	library(Rhpc)
	library(plyr)
	library(dplyr)
	library(ggplot2)
	library(reshape2)
	library(parallel)
	return(detectCores())
}

cores <- add.lib()

source("../indiff/indifference.r")
indiff <- round(indiff,2)
#indiff contains the indifference points of the 9 lotteries


# Grab our data
load("../data/agg-dat/Agg1Mil-S10k.Rda")
MM <- tbl_df(MM)

# Don't always want to use it all, there is tons of data
sam.prop <- 1

# Bounds for the means of data
lbound <- -1.9
ubound <- 1.55

# Certain data have rm inside the relevent range of the HL-MPL, others don't.
IN  <- MM %>%
	sample_frac(sam.prop) %>%
	filter(rm > lbound, rm < ubound) %>%
	filter(um < 1) %>%
	mutate(id = 1:n())

#OUT <- MM %>%
#	sample_frac(sam.prop) %>%
#	filter(rm < lbound | rm > ubound) %>%
#	mutate(id = 1:n())

# We no longer need the full dataset, it just wastes memory
rm(MM)

# Text naming the data.frame to use
to.USE <- "IN"

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

# Faceting, 
USE.w <- tbl_df(melt(get(to.USE), measure.vars=c("M.WP","V.WP","M.WC","V.WC") ))
USE.e <- tbl_df(melt(get(to.USE), measure.vars=c("M.EE","V.EE","M.WE0","V.WE0") ))

USE.w$variable  <-  ifelse(USE.w$variable == "M.WP" , "A) Mean Expected Welfare Proportion",
					ifelse(USE.w$variable == "V.WP" , "B) Var. of Expected Welfare Proportion",
					ifelse(USE.w$variable == "M.WC" , "C) Mean Expected Welfare Surplus",
													  "D) Var. of Expected Welfare Surplus") ))
USE.e$variable  <-  ifelse(USE.e$variable == "M.EE" , "A) Mean Expected Number of Errors",
					ifelse(USE.e$variable == "V.EE" , "B) Var. of Expected Number of Errors",
					ifelse(USE.e$variable == "M.WE0", "C) Mean Expected Proportion of No Error Choices",
													  "D) Var. of Expected Proportion of No Error Choices") ))

USE.w$variable <- factor(USE.w$variable, levels=c(	"A) Mean Expected Welfare Proportion",
													"B) Var. of Expected Welfare Proportion",
													"C) Mean Expected Welfare Surplus",
													"D) Var. of Expected Welfare Surplus"))
USE.e$variable <- factor(USE.e$variable, levels=c(	"A) Mean Expected Number of Errors",
													"B) Var. of Expected Number of Errors",
													"C) Mean Expected Proportion of No Error Choices",
													"D) Var. of Expected Proportion of No Error Choices"))

# Aggregate the two melded datasets to reduce RAM
USE.e <- USE.e %>%
	select(rm,rs,um,us,variable,value)
													
USE.w <- USE.w %>%
	select(rm,rs,um,us,variable,value)

USE <- USE.w %>%
	select(rm,rs,um)

USE$W.value <- USE.w$value
USE$E.value <- USE.e$value

USE$W.variable <- USE.w$variable
USE$E.variable <- USE.e$variable

m.den <- dnorm(USE$rm,mean=USE$rm,sd=USE$rs)

dists <- lapply(indiff,function(x){
	dnorm(x,mean=USE$rm,sd=USE$rs) / m.den
})
		
USE$den <- rowSums(do.call(cbind,dists))

# We now no longer need the "*USE*" dataframe, again it just wastes ram to keep it around
rm(list=c("USE.w","USE.e","to.USE"))

# What are the parameter names
names   <- c("den")
s.names <- c("um")

p.t <- c(rep("E",length(names)),rep("W",length(names)))

names   <- rep(names,2)
s.names <- rep(s.names,2)

to.plot <- data.frame(rbind(names,s.names,p.t))

getPlotted <- function(plot){

	add.lib()
    
    x.par <- as.character(plot[1])
    s.par <- as.character(plot[2])
    dtype <- as.character(plot[3])

	# Reference point for scale
	mid.ref <- (max(USE[[s.par]]) - ((max(USE[[s.par]])-(min(USE[[s.par]])))/2))
	limits <- c(min(USE[[s.par]]), max(USE[[s.par]]))

	q.dist <- (limits[2] - mid.ref) / 2

	bbr <- mid.ref - q.dist
	tbr <- mid.ref + q.dist

	breaks <- c(bbr,mid.ref,tbr,limits[2])
	breaks <- floor(breaks * 100) / 100
	ll <- ceiling(limits[1]*100)/100
	breaks <- c(ll,breaks)

	leg.title.size  <- 30
	leg.title.vjust <- -100

	leg.text.size   <- 24
	leg.text.angle  <- 22.5
	leg.position    <- "right"
	leg.direction	<- "horizontal"
	leg.direction	<- "vertical"
	leg.key.height  <- 1
	leg.key.width   <- .25

	x.title <- "Sum of relative densities"
			   
	x.title.size <- 36

	leg.title <- ifelse(s.par=="rm", "Mean of CRRA\n\n",
                 ifelse(s.par=="rs", "Standard Deviation\nof CRRA\n\n",
                 ifelse(s.par=="um", "Mean of Lambda\n\n","Standard Deviation\nof Lamda\n\n")))

	fac.size <- 36

	axis.text.size <- 24

	axis.x.angle <- 45
	axis.x.vjust <- .5

	axis.y.angle <- 0
	axis.y.vjust <- 0

	# Split the dataset up in a few ways to get multiple smoothed lines
	s.num <- 10

	USE.n <- USE

	USE.n <- USE.n %>% 
			mutate_(.dots=setNames(paste0("ntile(",s.par,",",s.num,")"),"bin")) %>%
			select_(.dots=list(x.par,s.par,"bin",paste0(dtype,".value"),paste0(dtype,".variable"))) 

	label <- c()
	
	for(i in 1:s.num){
		lower <- USE.n %>%
			filter(bin == i) %>%
			select_(.dots=list(s.par)) %>%
			min %>%
			round(2)
	
		upper <- USE.n %>%
			filter(bin == i) %>%
			select_(.dots=list(s.par)) %>%
			max %>%
			round(2)

		label <- c(label, paste0(lower," to ",upper))
	}

	# Split the dataset up in a few ways to get multiple smoothed lines
	USE.n$bin <- factor(USE.n$bin,levels=as.character(1:s.num),labels=label,ordered=T)

	p <- ggplot(data=USE.n,aes_string(x=x.par,y=paste0(dtype,".value"),color="bin"))

	p <- p + geom_point(alpha=.1)

	p <- p + scale_color_discrete(name=leg.title)  

	p <- p + guides(colour = guide_legend(override.aes = list(size = 10,alpha=1)))

	p <- p + facet_wrap( facets=paste0(dtype,".variable"), ncol=2, scale="free_y")

	p <- p + labs(x=x.title,y=NULL)

	p <- p + theme(legend.title=element_text(size=leg.title.size,vjust=.9),
					legend.text=element_text(size=leg.text.size,angle=leg.text.angle),
					legend.position=leg.position,
					legend.direction=leg.direction,
					legend.key.height=unit(leg.key.height,"in"),
					legend.key.width=unit(leg.key.width,"in"),
					axis.title.x=element_text(size=x.title.size),
					axis.text.y=element_text(size=axis.text.size,angle=axis.y.angle),
					axis.text.x=element_text(size=axis.text.size,angle=axis.x.angle,vjust=axis.x.vjust),
					strip.text=element_text(size=fac.size)
					)


	# Return the plot to the list
	dat <- ifelse(dtype=="W","Wel","Err")

	# Diminetions of plots
	h <- 8.50 - (.79*2) - .5	# Page is 8.5 x 11 inches, margins are .79 inches
	w <- 11.0 - (.79*2)

	save.scale <- 3

	fname <- paste("../data/agg-plots/D-",dat,"-points.jpg",sep="")

	print(paste(system("hostname"),"has started",fname))
	ggsave(filename=fname, plot=p, width=w, height=h, units="in",scale=save.scale )

	return(NULL)

}

# Plot generation is done over MPI, the cost of exporting objects is far outweighed
# by being able to run the calculations quicker

if(use.par){

	export(T,"add.lib","indiff","colors","USE","to.plot","getPlotted")

	c.lapply(X=to.plot,FUN=getPlotted)

	c.done()

} else{
	mclapply(to.plot,getPlotted,mc.cores=cores)
}


