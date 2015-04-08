# For ironing out bugs

rm(list=ls())
source("GBBfunctions.R")
source("GlobalParameters.R")


system.time(temp<-run.sim(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	initGens=initGens,
	Ve=Ve,
	mu=mu))



