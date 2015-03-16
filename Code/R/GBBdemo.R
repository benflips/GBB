rm(list=ls())
source("GBBfunctions.R")

# define global parameters
n<-100
nLoci<-10
maxGDisp<-4
eSize<-maxGDisp/nLoci
Ve<-(maxGDisp/4)^2


temp<-run.sim(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Ve=Ve)