# Need to find a way to aggregate data in a way that is sensible

# Grab a variable with all the filenames listed, this will be looped through
# later
files<-list.files(path=".", pattern="*.Rda", full.names=TRUE, recursive=TRUE)

# What aggregation do I want? Ultimately I want to see how welfare changes
# across distributions. So maybe multiply the probability of a pattern by the
# metrics I collected. Also, let's poke at the behavioral guys by also getting
# a metric on the proportion of expected subject that will violate
# deterministic EUT, and see how this changes across distributions.



PerFile <- function(file){

	load(file)

	# The dataset in the files were all saved as "D"

	W.EE <- D$PC * D$EE
	W.WP <- D$PC * D$WP
	W.WC <- D$PC * D$WC
	W.WE <- D$PC * D$WE

	W.NE <- D$PC * D$E.0 / 10000

	C <- apply(cbind(W.EE,W.WP,W.WC,W.WE,W.NE),2,sum)

	return(c(C,D$rmean[1],D$rstd[1],D$mumean[1],D$mustd[1]))

}

B <- apply(cbind(files),1,PerFile)

A <- data.frame(matrix(data=B, byrow=TRUE, ncol=nrow(B)))

colnames(A) <- c("W.EE","W.WP","W.WC","W.WE","W.NE","r.mean","r.std","mu.mean","mu.std")

# That was much easier than I thought it would be :)


	



