# Runs the basic spread model and collects populations
rm(list=ls())
source("GBBfunctions.R")
source("GlobalParameters.R")


parameters<-list(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	initGens=initGens)

spreadReps<-vector(mode="list", length=baseReps)

for (rr in 1:baseReps){
	cat("Replicate", rr, "...\n")
	spreadReps[[rr]]<-run.sim(n=n, 
		nLoci=nLoci, 
		eSize=eSize,
		Rmax=Rmax,
		Nstar=Nstar,
		k=k,
		initGens=initGens,
		Ve=Ve,
		mu=mu)
	save(parameters, spreadReps, file="../../Outputs/BaseSpread.RData")
}


