# Clear All
rm(list=ls())

library(ggplot2)
library(grid)
library(reshape2)
library(parallel)

cores <- detectCores()

## Haven't use, but will likely use
library(gridExtra)
library(np)

# Grab our data
load("../data/agg-dat/Agg1Mil-S10k.Rda")

# Don't always want to use it all, there is tons of data
sam <- runif(nrow(MM))
sam.prop <- .01
sam <- sam < sam.prop

M0 <- MM[which(sam),]

# Do I want to include some standard deviations above or below the bounds?
sd.include <- 0

# Functionally determine the alpha channel for the points
alph <- (1 / (exp(sam.prop))) 
alph <- .4

# Bounds for the means of data
lbound <- -1.7134
ubound <- 1.3684

lbound <- -1.9
ubound <- 1.55

# Certain data have rm inside the relevent range of the HL-MPL, others don't.
# Make the distiction
IN   <- M0[which(M0$rm > (lbound - sd.include*M0$rs) & M0$rm < (ubound + sd.include*M0$rs )),]
OUT  <- M0[which(M0$rm < (lbound + sd.include*M0$rs) | M0$rm > (ubound - sd.include*M0$rs )),]

IN$id  <- 1:nrow(IN)
OUT$id <- 1:nrow(OUT)

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
USE.m <- melt(get(USE), measure.vars=c("M.WP","M.WC","M.EE","M.WE0") )
USE.s <- melt(get(USE), measure.vars=c("V.WP","V.WC","V.EE","V.WE0") )

USE.m$variable  <-  ifelse(USE.m$variable == "M.WP" , "Mean Expected Welfare Proportion",
					ifelse(USE.m$variable == "M.WC" , "Mean Expected Welfare Surplus",
					ifelse(USE.m$variable == "M.EE" , "Mean Expected Errors",
					ifelse(USE.m$variable == "M.EE" , "Mean Expected Errors",
					ifelse(USE.m$variable == "M.WE0", "Mean Expected Zero Errors",USE.m$variable)
					))))

USE.s$variable  <-  ifelse(USE.s$variable == "V.WP" , "Variance of Expected Welfare Proportion",
					ifelse(USE.s$variable == "V.WC" , "Variance of Expected Welfare Surplus",
					ifelse(USE.s$variable == "V.EE" , "Variance of Expected Errors",
					ifelse(USE.s$variable == "V.EE" , "Variance of Expected Errors",
					ifelse(USE.s$variable == "V.WE0", "Variance of Expected Zero Errors",USE.s$variable)
					))))



# What are the parameter names
names <- c("rm","rs","um","us")
s.names <- c("rs","rm","us","um")

# Make a list to hold the plots for means
p.m <- list()
p.v <- list()

for(i in 1:length(names)){
    
	# Reference point for scale
	mid.ref <- (max(USE.m[[s.names[i]]]) - ((max(USE.m[[s.names[i]]])-(min(USE.m[[s.names[i]]])))/2))
	limits <- c(min(USE.m[[s.names[i]]]), max(USE.m[[s.names[i]]]))

	q.dist <- (limits[2] - mid.ref) / 2

	bbr <- mid.ref - q.dist
	tbr <- mid.ref + q.dist

	breaks <- c(bbr,mid.ref,tbr,limits[2])
	breaks <- floor(breaks * 100) / 100
	ll <- ceiling(limits[1]*100)/100
	breaks <- c(ll,breaks)


	# Values for the various plots

	leg.title.size  <- 22
	leg.title.vjust <- 0

	leg.text.size   <- 14
	leg.text.angle  <- 45
	leg.position    <- "top"
	leg.key.height  <- .25
	leg.key.width   <- 1

	x.title <- ifelse(names[i]=="rm", "Mean of CRRA",
			   ifelse(names[i]=="rs", "Standard Deviation of CRRA",
			   ifelse(names[i]=="um", "Mean of Lambda","Standard Deviation of Lamda")))

	x.title.size <- 24


	leg.title <- ifelse(s.names[i]=="rm", "Mean of CRRA",
                 ifelse(s.names[i]=="rs", "Standard Deviation\n of CRRA",
                 ifelse(s.names[i]=="um", "Mean of Lambda","Standard Deviation\n of Lamda")))


	fac.size <- 24


	# Plot the means
	p.m[[names[i]]] <- ggplot(data=USE.m, aes_string(x=names[i]))
	p.m[[names[i]]] <- p.m[[names[i]]] + scale_color_gradientn(name=leg.title,space="Lab",colours=colors,limits=limits,breaks=breaks) + 
										  scale_fill_gradientn(name=leg.title,space="Lab",colours=colors,limits=limits,breaks=breaks)
	p.m[[names[i]]] <- p.m[[names[i]]] + geom_point(shape=shape, alpha=alph, aes_string(y="value",color=s.names[i],fill=s.names[i]) )
	p.m[[names[i]]] <- p.m[[names[i]]] + facet_wrap( facets=~variable, scale="free_y")
	p.m[[names[i]]] <- p.m[[names[i]]] + labs(x=x.title,y=NULL)
	p.m[[names[i]]] <- p.m[[names[i]]] + theme(legend.title=element_text(size=leg.title.size,vjust=leg.title.vjust),
											   legend.text=element_text(size=leg.text.size,angle=leg.text.angle),
											   legend.position=leg.position,
											   legend.key.height=unit(leg.key.height,"in"),
											   legend.key.width=unit(leg.key.width,"in"),
											   axis.title.x=element_text(size=x.title.size),
											   strip.text=element_text(size=fac.size)
											   )

	# Plot the variances
	p.v[[names[i]]] <- ggplot(data=USE.s, aes_string(x=names[i]))
	p.v[[names[i]]] <- p.v[[names[i]]] + scale_color_gradientn(name=leg.title,space="Lab",colours=colors,limits=limits,breaks=breaks) + 
										  scale_fill_gradientn(name=leg.title,space="Lab",colours=colors,limits=limits,breaks=breaks)
	p.v[[names[i]]] <- p.v[[names[i]]] + geom_point(shape=shape, alpha=alph, aes_string(y="value",color=s.names[i],fill=s.names[i]) )
	p.v[[names[i]]] <- p.v[[names[i]]] + facet_wrap( facets=~variable, scale="free_y")
	p.v[[names[i]]] <- p.v[[names[i]]] + labs(x=x.title,y=NULL)
	p.v[[names[i]]] <- p.v[[names[i]]] + theme(legend.title=element_text(size=leg.title.size,vjust=leg.title.vjust),
											   legend.text=element_text(size=leg.text.size,angle=leg.text.angle),
											   legend.key.height=unit(leg.key.height,"in"),
											   legend.key.width=unit(leg.key.width,"in"),
											   legend.position=leg.position,
											   axis.title.x=element_text(size=x.title.size),
											   strip.text=element_text(size=fac.size)
											   )

}


# Save plots
saver <- function(pt) {

	p.t  <- pt[1]
	name <- pt[2]

	# Diminetions of plots
	w <- 6.92
	h <- (9.42/2) - .25

	save.scale <- 3

	fname <- paste("../data/agg-plots/",p.t,"-",name,".png",sep="")

	ggsave(filename=fname, plot=get(as.character(pt[1]))[[as.character(pt[2])]], width=w, height=h, units="in",scale=save.scale )

}

p.types <- c( rep("p.m",length(names)), rep("p.v",length(names)) )
p.types <- rbind(p.types , rep(names,2) )
p.types <- data.frame(p.types)

#rm(IN,OUT,M0,MM,USE.m,USE.s)

mclapply(p.types,saver,mc.cores=cores)
#lapply(p.types,saver)

