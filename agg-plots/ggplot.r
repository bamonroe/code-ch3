# Clear All
rm(list=ls())

library(ggplot2)
library(reshape2)
library(np)

load("../data/agg-dat/Agg1Mil-S10k.Rda")

sam <- runif(nrow(MM))
sam.prop <- .1
sam <- sam < sam.prop

M0 <- MM[which(sam),]

sd.include <- 0

alph <- (1 / (exp(sam.prop))) * .75

IN   <- M0[which(M0$rm > (-1.7134 - sd.include*M0$rs) & M0$rm < (1.3684 + sd.include*M0$rs )),]
OUT  <- M0[which(M0$rm < (-1.7134 + sd.include*M0$rs) | M0$rm > (1.3684 - sd.include*M0$rs )),]

# Colors
cy  <-"#32b0e6" # Cyan
yo  <-"#e6c832" # Yellow Orange
ma  <-"#e63271" # Magenta

# Text naming the data.frame to use
USE <- "IN"

# Reference point for scale
mid.ref <- (max(get(USE)$rs) - ((max(get(USE)$rs)-(min(get(USE)$rs)))/2))
limits <- c(min(get(USE)$rs), max(get(USE)$rs))
breaks <- c(limits[1],mid.ref,limits[2])

# Skeleton
p <- ggplot(data=get(USE), aes(x=rm))
# Change scale stuff
p <- p + scale_color_gradientn(space="Lab",colours=c(cy,yo,ma),limits=limits,breaks=breaks)
p <- p + scale_fill_gradientn(space="Lab",colours=c(cy,yo,ma),limits=limits,breaks=breaks)

# Collect the plots
rm.p <- list()

rm.p[["M.WP"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=M.WP) )
rm.p[["V.WP"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=V.WP) )
rm.p[["M.WC"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=M.WC) )
rm.p[["V.WC"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=V.WC) )
rm.p[["M.EE"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=M.EE) )
rm.p[["V.EE"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=V.EE) )
rm.p[["M.WE0"]] <- p + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=M.WE0) )
rm.p[["V.WE0"]] <- p + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=V.WE0) )

# Reference point for scale
mid.ref <- (max(get(USE)$rm) - ((max(get(USE)$rm)-(min(get(USE)$rm)))/2))
limits <- c(min(get(USE)$rm), max(get(USE)$rm))
breaks <- c(limits[1],mid.ref,limits[2])
# Skeleton
p <- ggplot(data=get(USE), aes(x=rs))
# Change scale stuff
p <- p + scale_color_gradientn(space="Lab",colours=c(cy,yo,ma),limits=limits,breaks=breaks)
p <- p + scale_fill_gradientn(space="Lab",colours=c(cy,yo,ma),limits=limits,breaks=breaks)

# Collect the plots
rs.p <- list()

rs.p[["M.WP"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rm,fill=rm,y=M.WP) )
rs.p[["V.WP"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rm,fill=rm,y=V.WP) )
rs.p[["M.WC"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rm,fill=rm,y=M.WC) )
rs.p[["V.WC"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rm,fill=rm,y=V.WC) )
rs.p[["M.EE"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rm,fill=rm,y=M.EE) )
rs.p[["V.EE"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=rm,fill=rm,y=V.EE) )
rs.p[["M.WE0"]] <- p + geom_point(shape=19, alpha=alph, aes(color=rm,fill=rm,y=M.WE0) )
rs.p[["V.WE0"]] <- p + geom_point(shape=19, alpha=alph, aes(color=rm,fill=rs,y=V.WE0) )

# Reference point for scale
mid.ref <- (max(get(USE)$us) - ((max(get(USE)$us)-(min(get(USE)$us)))/2))
limits <- c(min(get(USE)$us), max(get(USE)$us))
breaks <- c(limits[1],mid.ref,limits[2])
# Skeleton
p <- ggplot(data=get(USE), aes(x=um))
# Change scale stuff
p <- p + scale_color_gradientn(space="Lab",colours=c(cy,yo,ma),limits=limits,breaks=breaks)
p <- p + scale_fill_gradientn(space="Lab",colours=c(cy,yo,ma),limits=limits,breaks=breaks)
# Collect the plots
um.p <- list()

um.p[["M.WP"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=us,fill=us,y=M.WP) )
um.p[["V.WP"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=us,fill=us,y=V.WP) )
um.p[["M.WC"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=us,fill=us,y=M.WC) )
um.p[["V.WC"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=us,fill=us,y=V.WC) )
um.p[["M.EE"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=us,fill=us,y=M.EE) )
um.p[["V.EE"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=us,fill=us,y=V.EE) )
um.p[["M.WE0"]] <- p + geom_point(shape=19, alpha=alph, aes(color=us,fill=us,y=M.WE0) )
um.p[["V.WE0"]] <- p + geom_point(shape=19, alpha=alph, aes(color=us,fill=us,y=V.WE0) )

# Reference point for scale
mid.ref <- (max(get(USE)$um) - ((max(get(USE)$um)-(min(get(USE)$um)))/2))
limits <- c(min(get(USE)$um), max(get(USE)$um))
breaks <- c(limits[1],mid.ref,limits[2])
# Skeleton
p <- ggplot(data=get(USE), aes(x=us))
# Change scale stuff
p <- p + scale_color_gradientn(space="Lab",colours=c(cy,yo,ma),limits=limits,breaks=breaks)
p <- p + scale_fill_gradientn(space="Lab",colours=c(cy,yo,ma),limits=limits,breaks=breaks)
# Collect the plots
us.p <- list()

us.p[["M.WP"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=um,fill=um,y=M.WP) )
us.p[["V.WP"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=um,fill=um,y=V.WP) )
us.p[["M.WC"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=um,fill=um,y=M.WC) )
us.p[["V.WC"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=um,fill=um,y=V.WC) )
us.p[["M.EE"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=um,fill=um,y=M.EE) )
us.p[["V.EE"]]  <- p + geom_point(shape=19, alpha=alph, aes(color=um,fill=um,y=V.EE) )
us.p[["M.WE0"]] <- p + geom_point(shape=19, alpha=alph, aes(color=um,fill=um,y=M.WE0) )
us.p[["V.WE0"]] <- p + geom_point(shape=19, alpha=alph, aes(color=um,fill=um,y=V.WE0) )

