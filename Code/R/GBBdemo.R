#To do: run sims on current code to confirm expectations

rm(list=ls())
source("GBBfunctions.R")
source("GlobalParameters.R")



temp<-run.sim(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar)

plot(temp$N, type="l")
plot(temp$mean.di, type="l")



