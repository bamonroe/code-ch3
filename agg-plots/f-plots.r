# Clear All
rm(list=ls())

library(ggplot2)
library(reshape2)

## Haven't use, but will likely use
library(gridExtra)
library(np)


# Grab our data
load("../data/agg-dat/Agg1Mil-S10k.Rda")

# Don't always want to use it all, there is tons of data
sam <- runif(nrow(MM))
sam.prop <- 1
sam <- sam < sam.prop

M0 <- MM[which(sam),]

# Do I want to include some standard deviations above or below the bounds?
sd.include <- 0

# Functionally determine the alpha channel for the points
alph <- (1 / (exp(sam.prop))) 
alph <- .1

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
shape <- 19

# Faceting, 
USE.m <- melt(get(USE), measure.vars=c("M.WP","M.WC","M.EE","M.WE0") )
USE.s <- melt(get(USE), measure.vars=c("V.WP","V.WC","V.EE","V.WE0") )

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

	breaks <- c(limits[1],bbr,mid.ref,tbr,limits[2])

	# Plot the means
	p.m[[names[i]]] <- ggplot(data=USE.m, aes_string(x=names[i]))
	p.m[[names[i]]] <- p.m[[names[i]]] + scale_color_gradientn(space="Lab",colours=colors,limits=limits,breaks=breaks) + scale_fill_gradientn(space="Lab",colours=colors,limits=limits,breaks=breaks)
	p.m[[names[i]]] <- p.m[[names[i]]] + geom_point(shape=shape, alpha=alph, aes_string(y="value",color=s.names[i],fill=s.names[i]) )
	p.m[[names[i]]] <- p.m[[names[i]]] + facet_wrap( facets=~variable, scale="free_y")

	# Plot the variances
	p.v[[names[i]]] <- ggplot(data=USE.s, aes_string(x=names[i]))
	p.v[[names[i]]] <- p.v[[names[i]]] + scale_color_gradientn(space="Lab",colours=colors,limits=limits,breaks=breaks) + scale_fill_gradientn(space="Lab",colours=colors,limits=limits,breaks=breaks)
	p.v[[names[i]]] <- p.v[[names[i]]] + geom_point(shape=shape, alpha=alph, aes_string(y="value",color=s.names[i],fill=s.names[i]) )
	p.v[[names[i]]] <- p.v[[names[i]]] + facet_wrap( facets=~variable, scale="free_y")

}

