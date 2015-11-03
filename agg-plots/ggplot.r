# Clear All
rm(list=ls())


library(ggplot2)
library(reshape2)
library(np)


#load("../data/agg-dat/Agg30k-S500.Rda")
load("../data/agg-dat/Agg1Mil-S10k.Rda")

sam <- runif(nrow(MM))
sam <- sam < 1

M0 <- MM[which(sam),]

sd.include <- 0


IN   <- M0[which(M0$rm > (-1.7134 - sd.include*M0$rs) & M0$rm < (1.3684 + sd.include*M0$rs )),]
OUT  <- M0[which(M0$rm < (-1.7134 + sd.include*M0$rs) | M0$rm > (1.3684 - sd.include*M0$rs )),]

rm <- ggplot(data=IN, aes(x=rm))
rm.p <- list()

rm.p[["M.WP"]]  <- rm + geom_point(shape=19, color="black", alpha=.2, aes(y=M.WP) )
rm.p[["V.WP"]]  <- rm + geom_point(shape=19, color="black", alpha=.2, aes(y=V.WP) )
rm.p[["M.WC"]]  <- rm + geom_point(shape=19, color="black", alpha=.2, aes(y=M.WC) )
rm.p[["V.WC"]]  <- rm + geom_point(shape=19, color="black", alpha=.2, aes(y=V.WC) )
rm.p[["M.EE"]]  <- rm + geom_point(shape=19, color="black", alpha=.2, aes(y=M.EE) )
rm.p[["V.EE"]]  <- rm + geom_point(shape=19, color="black", alpha=.2, aes(y=V.EE) )
rm.p[["M.WE0"]] <- rm + geom_point(shape=19, color="black", alpha=.2, aes(y=M.WE0) )
rm.p[["V.WE0"]] <- rm + geom_point(shape=19, color="black", alpha=.2, aes(y=V.WE0) )


rm.p[["M.WP"]]  <- rm + geom_point(shape=19, alpha=.2, aes(color=rs,y=M.WP) )
rm.p[["V.WP"]]  <- rm + geom_point(shape=19, alpha=.2, aes(color=rs,y=V.WP) )
rm.p[["M.WC"]]  <- rm + geom_point(shape=19, alpha=.2, aes(color=rs,y=M.WC) )
rm.p[["V.WC"]]  <- rm + geom_point(shape=19, alpha=.2, aes(color=rs,y=V.WC) )
rm.p[["M.EE"]]  <- rm + geom_point(shape=19, alpha=.2, aes(color=rs,y=M.EE) )
rm.p[["V.EE"]]  <- rm + geom_point(shape=19, alpha=.2, aes(color=rs,y=V.EE) )
rm.p[["M.WE0"]] <- rm + geom_point(shape=19, alpha=.2, aes(color=rs,y=M.WE0) )
rm.p[["V.WE0"]] <- rm + geom_point(shape=19, alpha=.2, aes(color=rs,y=V.WE0) )


rs <- ggplot(data=IN, aes(x=rs))
rs.p <- list()

rs.p[["M.WP"]]  <- rs + geom_point(shape=19,  alpha=.2, aes(color=rm,y=M.WP) )
rs.p[["V.WP"]]  <- rs + geom_point(shape=19,  alpha=.2, aes(color=rm,y=V.WP) )
rs.p[["M.WC"]]  <- rs + geom_point(shape=19,  alpha=.2, aes(color=rm,y=M.WC) )
rs.p[["V.WC"]]  <- rs + geom_point(shape=19,  alpha=.2, aes(color=rm,y=V.WC) )
rs.p[["M.EE"]]  <- rs + geom_point(shape=19,  alpha=.2, aes(color=rm,y=M.EE) )
rs.p[["V.EE"]]  <- rs + geom_point(shape=19,  alpha=.2, aes(color=rm,y=V.EE) )
rs.p[["M.WE0"]] <- rs + geom_point(shape=19,  alpha=.2, aes(color=rm,y=M.WE0) )
rs.p[["V.WE0"]] <- rs + geom_point(shape=19,  alpha=.2, aes(color=rm,y=V.WE0) )


um <- ggplot(data=IN, aes(x=um))
um.p <- list()

um.p[["M.WP"]]  <- um + geom_point(shape=19, color="black", alpha=.2, aes(y=M.WP) )
um.p[["V.WP"]]  <- um + geom_point(shape=19, color="black", alpha=.2, aes(y=V.WP) )
um.p[["M.WC"]]  <- um + geom_point(shape=19, color="black", alpha=.2, aes(y=M.WC) )
um.p[["V.WC"]]  <- um + geom_point(shape=19, color="black", alpha=.2, aes(y=V.WC) )
um.p[["M.EE"]]  <- um + geom_point(shape=19, color="black", alpha=.2, aes(y=M.EE) )
um.p[["V.EE"]]  <- um + geom_point(shape=19, color="black", alpha=.2, aes(y=V.EE) )
um.p[["M.WE0"]] <- um + geom_point(shape=19, color="black", alpha=.2, aes(y=M.WE0) )
um.p[["V.WE0"]] <- um + geom_point(shape=19, color="black", alpha=.2, aes(y=V.WE0) )


us <- ggplot(data=IN, aes(x=us))
us.p <- list()

us.p[["M.WP"]]  <- us + geom_point(shape=19, color="black", alpha=.2, aes(y=M.WP) )
us.p[["V.WP"]]  <- us + geom_point(shape=19, color="black", alpha=.2, aes(y=V.WP) )
us.p[["M.WC"]]  <- us + geom_point(shape=19, color="black", alpha=.2, aes(y=M.WC) )
us.p[["V.WC"]]  <- us + geom_point(shape=19, color="black", alpha=.2, aes(y=V.WC) )
us.p[["M.EE"]]  <- us + geom_point(shape=19, color="black", alpha=.2, aes(y=M.EE) )
us.p[["V.EE"]]  <- us + geom_point(shape=19, color="black", alpha=.2, aes(y=V.EE) )
us.p[["M.WE0"]] <- us + geom_point(shape=19, color="black", alpha=.2, aes(y=M.WE0) )
us.p[["V.WE0"]] <- us + geom_point(shape=19, color="black", alpha=.2, aes(y=V.WE0) )



