path <- "../data/agg-dat/"

fnames <- dir(path,pattern=".Rda")

for(f in 1: length(fnames)){
    
    load(paste(path,fnames[f],sep=""))

	P <- data.frame(rm=MM$rm,rs=MM$rs,um=MM$um,us=MM$us)

	rm(MM)

	save(P,file=paste(path,"P-",fnames[f],sep=""))

}
