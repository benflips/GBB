# examines the effectiveness of constant barrier to varying introduction extents.
rm(list=ls())
source("GBBfunctions.R")
source("GlobalParameters.R")

ext.parameters<-list(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	initGens=initGens,
	Ve=Ve, 
	mu=mu, 
	defBar=defBar, 
	lead=lead,
	bbMonitorGens=bbMonitorGens)
	
extSeq<-seq(0, lead, 2)
extentMat<-cbind(rep(extSeq, times=rep(bbReps, length(extSeq))), breachTime=rep(NA, length(extSeq)*bbReps))


for (ee in 1:length(extSeq)){
	cat("Extent", extSeq[ee], "\n")
	for (rr in 1:bbReps){
		cat("\tReplicate", rr, "\n")
		extentMat[((ee-1)*bbReps+rr),2]<-run.sim.barr(n, 
			nLoci, 
			eSize, 
			Rmax, 
			Nstar, 
			k, 
			initGens, 
			Ve, 
			mu, 
			defBar, 
			extent=extSeq[ee], 
			lead,
			bbMonitorGens)
		save(ext.parameters, extentMat, file="../../Outputs/ExtentTests.RData")
	}
	if (sum(extentMat[extentMat[,1]==extSeq[ee],2]==bbMonitorGens)==bbReps) break
}

#save(ext.parameters, extentMat, file="../../Outputs/ExtentTests.RData")