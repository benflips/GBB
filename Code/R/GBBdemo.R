# For ironing out bugs
# and profiling for speed.
rm(list=ls())

library(aprof)


source("GBBfunctions.R")
source("GlobalParameters.R")

initGens<-20

tmp<-tempfile()
Rprof(tmp,line.profiling=TRUE)

temp<-run.sim(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	initGens=initGens,
	Ve=Ve,
	mu=mu)

Rprof(append=FALSE)

fooaprof<-aprof("GBBfunctions.R",tmp)
#display basic information, summarize and plot the object

	 fooaprof
	 summary(fooaprof)
	 plot(fooaprof)
	 profileplot(fooaprof)


