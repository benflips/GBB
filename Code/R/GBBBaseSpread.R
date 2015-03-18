rm(list=ls())
source("GBBfunctions.R")
source("GlobalParameters.R")


spreadReps<-vector(mode="list", length=baseReps)

for (rr in 1:baseReps){
	spreadReps[[rr]]<-run.sim(n=n, 
		nLoci=nLoci, 
		eSize=eSize,
		Rmax=Rmax,
		Nstar=Nstar,
		k=k,
		initGens=initGens)
}
