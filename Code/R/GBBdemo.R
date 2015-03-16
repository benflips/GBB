rm(list=ls())
source("GBBfunctions.R")

# define global parameters
n<-100
nLoci<-10
Rmax<-5
Nstar<-30
maxGDisp<-log(4)
eSize<-maxGDisp/nLoci



temp<-run.sim(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax)
plot(temp)

hist(temp$x)
