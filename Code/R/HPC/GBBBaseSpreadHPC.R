source.dir<-"/home/jc227089/GenBackBurn/GBB/Code/R/HPC/"
outdir<-"/home/jc227089/GenBackBurn/Outputs"
source(paste(sourcedir, "GBBfunctions.R", sep=""))
source(paste(sourcedir, "GlobalParameters.R", sep=""))

args=(commandArgs(TRUE))

#evaluate the arguments
# input argument will repNumber
for(i in 1:length(args)) {
	 eval(parse(text=args[[i]]))
}



spreadRep<-run.sim(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	initGens=initGens)

parameters<-list(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	initGens=initGens)

save(parameters, spreadReps, 
	file=paste(outdir, "BaseSpread", repNumber, ".RData", sep=""))
