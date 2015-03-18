# For ironing out bugs

rm(list=ls())
source("GBBfunctions.R")
source("GlobalParameters.R")


initGens<-40

temp<-run.sim(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	initGens=initGens)

plot(temp$N, type="l")
plot(temp$mean.di, type="l")
plot(temp$X.g, type="l")
plot(temp$pop[,"X"], temp$pop[,"di"])

plot(tapply(temp$pop[,"di"], (temp$pop[,"X"] %/% 1), mean), type="l")
plot(temp$pop[,"X"], dens(temp$pop))


