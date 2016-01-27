library(Rhpc)

to.export <- c()

export <- function(obj){

	if(is.logical(obj)){
	
		Rhpc_export(cl,to.export)
	}else{
		if( length(to.export) == 0 ) {
			to.export <<- c(obj)
		}else{
			to.export <<- c(to.export,obj)
		}
	}
}

m.lapply <- function(X,FUN){
	Rhpc_lapply(cl,X=X,FUN=FUN)
}

# Rhcp requires the total worker number, not the total CPU number
worker.num <- read("/borg/cpu.num") - 1

Rhpc_initialize()
cl <- Rhpc_getHandle(worker.num)
