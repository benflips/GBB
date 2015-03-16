#To do: run sims on current code to confirm expectations

rm(list=ls())
source("GBBfunctions.R")

# define global parameters
n<-10
nLoci<-10
Rmax<-5
Nstar<-100
maxGDisp<-log(4)
eSize<-maxGDisp/nLoci



temp<-run.sim(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar)
plot(temp)
plot(temp$y, temp$z)
hist(temp$x)
