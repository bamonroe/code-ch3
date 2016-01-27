library(Rhpc)

to.export <- list()

export <- function(...){

	obj <- list(...)

	if(is.logical(obj[[1]])){
		# If the first element is a logical T then add all following elements
		# then export the to.export list
		if (obj[[1]] == T){

			if(length(obj)>1){
				for(i in 2:length(obj)){
					te.len <- length(to.export) + 1
					to.export[[te.len]] <<- obj[[i]]
				}
			} 

			for(i in 1:length(to.export)){
			
				# Check for character vectors and existance
				is.char <- ifelse(is.character(to.export[[i]]),T,F)
				it.exists <- ifelse( is.char == T && exists(to.export[[i]]),T,F )

				if( !is.char ) 	 stop(paste("Argument",i,"is not a string"))
				if( !it.exists ) stop(paste("Argument",i,"does not refer to an existing object"))

			}

			to.export <- do.call(c,to.export)

			# IF the value True is passed, conduct the export
			s.msg <- paste("Starting export of:",to.export)
			e.msg <- paste("Exported:",to.export)

			print("Starting export of:")
			print(to.export)
			Rhpc_Export(cluster.object,to.export)
			print("Exported:")
			return(print(to.export))

		} else if(obj[[1]] == F){
			# If the first element is a logical F, then clear the to.export
			# list and add the following elements back to it

			to.export <- list()
			if(length(obj)>1){
				for(i in 2:length(obj)){
					te.len <- length(to.export) + 1
					to.export[[te.len]] <<- obj[[i]]
				}
			} 
		}
	}

	for(i in 1:length(obj)){
	
		# Check for character vectors and existance
		is.char <- ifelse(is.character(obj[i]),T,F)
		it.exists <- ifelse( is.char == T && exists(obj[i]),T,F )

		if( !is.char ) 	 warning(paste("Argument",i,"is not a string"))
		if( !it.exists ) warning(paste("Argument",i,"does not refer to an existing object"))
	
		te.len <- length(to.export) + 1
		to.export[[te.len]] <<- obj[[i]]

	}

}

# Some pure wrappers
c.lapply <- function(X,FUN){
	Rhpc_lapply(cluster.object,X=X,FUN=FUN)
}
c.done <- function(){
	Rhpc_finalize()
}

# Rhcp requires the total worker number, not the total CPU number
worker.num <- as.integer(readChar("/borg/cpu.num",file.info("/borg/cpu.num")$size - 1)) - 1
# Initialize this cluster
Rhpc_initialize()
cluster.object <- Rhpc_getHandle(worker.num)
