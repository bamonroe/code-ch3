# Clear All
rm(list=ls())


library(ggplot2)
library(reshape2)
library(np)

load("../data/agg-dat/Agg1Mil-S10k.Rda")

sam <- runif(nrow(MM))
sam.prop <- 1
sam <- sam < sam.prop

M0 <- MM[which(sam),]

sd.include <- 0

alph <- (1 / (exp(sam.prop))) * .75

IN   <- M0[which(M0$rm > (-1.7134 - sd.include*M0$rs) & M0$rm < (1.3684 + sd.include*M0$rs )),]
OUT  <- M0[which(M0$rm < (-1.7134 + sd.include*M0$rs) | M0$rm > (1.3684 - sd.include*M0$rs )),]


# Colors
c  <-"#32b0e6" # Cyan
yo <-"#e6c832" # Yellow Orange
m  <-"#e63271" # Magenta


# Reference point for scale
mid <- (max(M0$rs) - ((max(M0$rs)-(min(M0$rs)))/2))
# Skeleton
rm <- ggplot(data=IN, aes(x=rm))
# Change scale stuff
#rm <- rm + scale_color_gradient2(low=c,mid=yo,high=m,midpoint= mid, )

rm <- rm + scale_color_gradient2(low=c,mid=yo,high=m,midpoint=mid)
rm <- rm + scale_fill_gradient2(low=c,mid=yo,high=m,midpoint=mid)

# Collect the plots
rm.p <- list()

rm.p[["M.WP"]]  <- rm + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=M.WP) )
rm.p[["V.WP"]]  <- rm + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=V.WP) )
rm.p[["M.WC"]]  <- rm + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=M.WC) )
rm.p[["V.WC"]]  <- rm + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=V.WC) )
rm.p[["M.EE"]]  <- rm + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=M.EE) )
rm.p[["V.EE"]]  <- rm + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=V.EE) )
rm.p[["M.WE0"]] <- rm + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=M.WE0) )
rm.p[["V.WE0"]] <- rm + geom_point(shape=19, alpha=alph, aes(color=rs,fill=rs,y=V.WE0) )

stop("here")

# Reference point for scale
mid <- (max(M0$rm) - ((max(M0$rm)-(min(M0$rm)))/2))
# Skeleton
rs <- ggplot(data=IN, aes(x=rs))
# Change scale stuff
rs <- rs + scale_color_gradient2(low=c,mid=yo,high=m,midpoint= mid )

# Collect the plots
rs.p <- list()

rs.p[["M.WP"]]  <- rs + geom_point(shape=19,  alpha=alph, aes(color=rm,y=M.WP) )
rs.p[["V.WP"]]  <- rs + geom_point(shape=19,  alpha=alph, aes(color=rm,y=V.WP) )
rs.p[["M.WC"]]  <- rs + geom_point(shape=19,  alpha=alph, aes(color=rm,y=M.WC) )
rs.p[["V.WC"]]  <- rs + geom_point(shape=19,  alpha=alph, aes(color=rm,y=V.WC) )
rs.p[["M.EE"]]  <- rs + geom_point(shape=19,  alpha=alph, aes(color=rm,y=M.EE) )
rs.p[["V.EE"]]  <- rs + geom_point(shape=19,  alpha=alph, aes(color=rm,y=V.EE) )
rs.p[["M.WE0"]] <- rs + geom_point(shape=19,  alpha=alph, aes(color=rm,y=M.WE0) )
rs.p[["V.WE0"]] <- rs + geom_point(shape=19,  alpha=alph, aes(color=rm,y=V.WE0) )


um <- ggplot(data=IN, aes(x=um))
um.p <- list()

um.p[["M.WP"]]  <- um + geom_point(shape=19, color="black", alpha=alph, aes(y=M.WP) )
um.p[["V.WP"]]  <- um + geom_point(shape=19, color="black", alpha=alph, aes(y=V.WP) )
um.p[["M.WC"]]  <- um + geom_point(shape=19, color="black", alpha=alph, aes(y=M.WC) )
um.p[["V.WC"]]  <- um + geom_point(shape=19, color="black", alpha=alph, aes(y=V.WC) )
um.p[["M.EE"]]  <- um + geom_point(shape=19, color="black", alpha=alph, aes(y=M.EE) )
um.p[["V.EE"]]  <- um + geom_point(shape=19, color="black", alpha=alph, aes(y=V.EE) )
um.p[["M.WE0"]] <- um + geom_point(shape=19, color="black", alpha=alph, aes(y=M.WE0) )
um.p[["V.WE0"]] <- um + geom_point(shape=19, color="black", alpha=alph, aes(y=V.WE0) )


us <- ggplot(data=IN, aes(x=us))
us.p <- list()

us.p[["M.WP"]]  <- us + geom_point(shape=19, color="black", alpha=alph, aes(y=M.WP) )
us.p[["V.WP"]]  <- us + geom_point(shape=19, color="black", alpha=alph, aes(y=V.WP) )
us.p[["M.WC"]]  <- us + geom_point(shape=19, color="black", alpha=alph, aes(y=M.WC) )
us.p[["V.WC"]]  <- us + geom_point(shape=19, color="black", alpha=alph, aes(y=V.WC) )
us.p[["M.EE"]]  <- us + geom_point(shape=19, color="black", alpha=alph, aes(y=M.EE) )
us.p[["V.EE"]]  <- us + geom_point(shape=19, color="black", alpha=alph, aes(y=V.EE) )
us.p[["M.WE0"]] <- us + geom_point(shape=19, color="black", alpha=alph, aes(y=M.WE0) )
us.p[["V.WE0"]] <- us + geom_point(shape=19, color="black", alpha=alph, aes(y=V.WE0) )



