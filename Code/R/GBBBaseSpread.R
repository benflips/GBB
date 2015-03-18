rm(list=ls())
source("GBBfunctions.R")
source("GlobalParameters.R")


spreadReps<-vector(mode="list", length=baseReps)

for (rr in 1:baseReps){
	cat("Replicate", rr, "...\n")
	spreadReps[[rr]]<-run.sim(n=n, 
		nLoci=nLoci, 
		eSize=eSize,
		Rmax=Rmax,
		Nstar=Nstar,
		k=k,
		initGens=initGens)
}

parameters<-list(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	initGens=initGens)

save(parameters, spreadReps, file="../../Outputs/BaseSpread.RData")
