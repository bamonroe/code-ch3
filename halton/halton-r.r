# Make a Halton Sequence for use with MSL

#Todo : The ananymous function is very simple. It may be worth writing in Rcpp

halton <- function(init,H,prime){
	sapply(c(init:H),function(index,base=prime){
    
		r <- 0
		f <- 1
		i <- index

		while(i > 0){
			f <- f / base
			r <- r + f*(i %% base)
			i <- floor(i /base)
		}
		return(r)
	})
}

