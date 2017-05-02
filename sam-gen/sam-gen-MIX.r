# Recreate Error Calculation
# Clear all
rm(list = ls())

#Libraries
library(microbenchmark)
library(Rcpp)
library(halton)
library(dplyr)

# C++ functions for this script
sourceCpp("../Rcpp/sam-gen.cpp")

# Pattern matrix, where each column is a choice pattern
patmat <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
patmat <- t(patmat)

# Set up our instrument
A0 <- rep(1.60,10)
A1 <- rep(2.00,10)
B0 <- rep(0.10,10)
B1 <- rep(3.85,10)

pA0 <- seq(.9,0,by=-.1)
pB0 <- seq(.9,0,by=-.1)
pA1 <- seq(.1,1,by=.1)
pB1 <- seq(.1,1,by=.1)

Max <- rep(3.85,10)
Min <- c(rep(0.10,9),2)

# How many simulations to run per sample
snum <- 500000

# How many choices per subject
cnum <- length(A0)

# Set up the Halton sequences
burn <- 30

# Set up the Halton sequences
Hr <- halton(snum, prime = 3)	# for r
Hu <- halton(snum, prime = 7)	# for mu

# Define the distributional parameters for EUT
rm <- .65
rs <- .3

um <- .35
us <- .4

# Fechner will use a gamma distribution, so need to back out shape and
# scale parameters
ku <- (um^2)/(us^2)
tu <- (us^2)/um

sim <- matrix( c(qnorm(Hr, mean = rm, sd = rs), 
								 qgamma(Hu, shape = ku, scale = tu)), 
								 byrow = TRUE, ncol = snum)

# Give the very wide bounds we've made for ourselves, occasionally we get 
# a fechner value equal to 0. This creates an "oh shi-" problem later on.
# Keep replacing 0s with another random draw from the same distribution until
# they go away.
#while(any(sim[,2]==0)){
#	sim[,2] <- ifelse(sim[,2]==0 , rgamma(1,shape = k, scale = t),sim[,2] )
#}
#while(anyNA(sim[,2])){
#	sim[,2] <- ifelse(is.na(sim[,2]) , rgamma(1,shape = k, scale = t),sim[,2] )
#}

# getRes is written in C++ and sourced from sam-gen.cpp
Res <- getRes(sim,A0,A1,B0,B1,
				pA0,pA1,pB0,pB1,
				Max,Min)


# Each column of Res is created from the following combined vector
#c(A.err , B.err , CEA , CEB , CEM , pA , pB)

M <- rowMeans(Res)

Errors <- Res[1:20, ]
Cert   <- Res[21:40, ]
CEMax  <- rbind(Res[41:50, ] , Res[41:50, ])
Probs  <- Res[51:70, ]

# DDcpp is written in C++ and sourced from sam-gen.cpp
EUT <- data.frame( DDcpp(pat=patmat,M=M,Errors=Errors,Cert=Cert,CEMax=CEMax,Probs=Probs) )

# Free up some memory
rm(Errors,Cert,CEMax,Probs)

# I choose to do this step in R to ensure precision. Highly vectorized
# operations are not very time consuming in R
EUT[,5:15] <- EUT[,5:15] / snum
EUT[,3:4]  <- EUT[,3:4] / cnum

# Add the pattern and parameters to D for clarity
EUT$pattern <- apply(patmat,2,toString)
EUT$rmean  <- rm
EUT$mumean <- um
EUT$rstd   <- rs
EUT$mustd  <- us

# Add column labels
colnames(EUT) <- c("EE","PC","WC","WP","E.0", "E.1", "E.2", "E.3", "E.4", "E.5", "E.6", "E.7", "E.8", "E.9", "E.10","Pattern","rm","rs","um","us")

# Sort by pattern, this is necessary for mixing 
EUT <- EUT[order(EUT$PC, decreasing=T),]

# Add row labels
rownames(EUT) <- 1:nrow(EUT)

# Sort by pattern for mixing
EUT <- EUT[order(EUT$Pattern, decreasing=T),]

## Now the RDU stuff
Ha <- halton(snum, prime = 11)	# for r
Hb <- halton(snum, prime = 13)	# for mu

# Fechner will use a gamma distribution, so need to back out shape and
# scale parameters

rm <- .65
rs <- .3

um <- .15
us <- .2

am <- 1.3
as <- .1

bm <- 1.2
bs <- .1

ku <- (um^2)/(us^2)
tu <- (us^2)/um

ka <- (am^2)/(as^2)
ta <- (as^2)/am

kb <- (bm^2)/(bs^2)
tb <- (bs^2)/bm

sim <- matrix( c(qnorm(Hr, mean = rm, sd = rs), 
								 qgamma(Hu, shape = ku, scale = tu),
								 qgamma(Ha, shape = ka, scale = ta),
								 qgamma(Hb, shape = kb, scale = tb)), 
								 byrow = TRUE, ncol = snum)

# getRes is written in C++ and sourced from sam-gen.cpp
Res <- getRes(sim,A0,A1,B0,B1,
				pA0,pA1,pB0,pB1,
				Max,Min)

# Each column of Res is created from the following combined vector
#c(A.err , B.err , CEA , CEB , CEM , pA , pB)

M <- rowMeans(Res)

Errors <- Res[1:20, ]
Cert   <- Res[21:40, ]
CEMax  <- rbind(Res[41:50, ] , Res[41:50, ])
Probs  <- Res[51:70, ]

# Free up some memory
rm(Res)

# DDcpp is written in C++ and sourced from sam-gen.cpp
RDU <- data.frame( DDcpp(pat=patmat,M=M,Errors=Errors,Cert=Cert,CEMax=CEMax,Probs=Probs) )

# Free up some memory
rm(Errors,Cert,CEMax,Probs)

# I choose to do this step in R to ensure precision. Highly vectorized
# operations are not very time consuming in R
RDU[,5:15] <- RDU[,5:15] / snum
RDU[,3:4] <- RDU[,3:4] / cnum

# Add the pattern and parameters to D for clarity
RDU$pattern <- apply(patmat,2,toString)
RDU$rmean  <- rm
RDU$mumean <- um
RDU$rstd   <- rs
RDU$mustd  <- us

# Add column labels
colnames(RDU) <- c("EE","PC","WC","WP","E.0", "E.1", "E.2", "E.3", "E.4", "E.5", "E.6", "E.7", "E.8", "E.9", "E.10","Pattern","rm","rs","um","us")

# Sort by probability, this makes the most sense
RDU <- RDU[order(RDU$PC, decreasing=T),]

# Add row labels
rownames(RDU) <- 1:nrow(RDU)

# Sort by pattern for mixing
RDU <- RDU[order(RDU$Pattern, decreasing=T),]


# Start Mixing
EUT.prop <- 0.7

MIX <- EUT[,1:15]*EUT.prop + RDU[,1:15]*(1-EUT.prop)
MIX$Pattern <- EUT$Pattern

MAG <- data.frame( EUTprop = c((EUT$PC*EUT.prop) / ((RDU$PC*(1-EUT.prop)) + (EUT$PC*EUT.prop))), 
									Pattern = MIX$Pattern)

# Sort Back to probabilities
EUT <- EUT[order(EUT$PC, decreasing=T),]
RDU <- RDU[order(RDU$PC, decreasing=T),]

MIX <- MIX[order(MIX$PC, decreasing=T),]
MIX$rank <- 1:nrow(MIX)
MIX <- MIX[order(MIX$Pattern, decreasing=T),]
MAG$rank <- MIX$rank

MIX <- MIX[order(MIX$PC, decreasing=T),]
MAG <- MAG[order(MAG$rank, decreasing=F),]
MIX$EUTprop <- MAG$EUTprop


EUT <- tbl_df(EUT)
RDU <- tbl_df(RDU)
MIX <- tbl_df(MIX)

print(head(EUT,n=10))
print(head(RDU,n=10))
print(head(MIX,n=10))

print(EUT.prop)

# Now lets get this all in a CSV that can be directly imported into latex

TopTenEUT <- EUT %>%
				select(Pattern, PC, EE, WP, WC, E.0, E.1) %>%
				top_n(10 , PC)


TopTenRDU <- RDU %>%
				select(Pattern, PC, EE, WP, WC, E.0, E.1) %>%
				top_n(10 , PC)

TopTenMIX <- MIX %>%
				select(Pattern, PC, EE, WP, WC, E.0, E.1) %>%
				top_n(10 , PC)

subpop <- c(1:10, "PC" , "EE", "WP", "WC", "E.0", "E.1")

for(mod in c("TopTenEUT", "TopTenRDU")){
	pat <- lapply( 1:10, function(i){
					get(mod)$Pattern[i]				%>%
					strsplit( split = c(",")) %>%
					do.call(what = c)					%>%
					trimws()									%>%
					as.numeric()
				})
	pat <- do.call(rbind, pat)

	out <- get(mod) %>%
						select(-Pattern)

	out <- cbind(pat,out)

	assign(paste0(mod,".out") , out)

	write.csv(out,file = paste0("../tables/",mod,".csv"), row.names = FALSE, quote = FALSE)

}

pat <- lapply( 1:10, function(i){
				TopTenMIX$Pattern[i]				%>%
				strsplit( split = c(",")) %>%
				do.call(what = c)					%>%
				trimws()									%>%
				as.numeric()
			})
pat <- do.call(rbind, pat)

out <- TopTenMIX %>%
					select(-Pattern, -E.1)

out <- cbind(pat, EUT.prop, out)

assign(paste0(mod,".out") , out)

write.csv(out,file = paste0("../tables/",mod,".csv"), row.names = FALSE, quote = FALSE)

