# examines the effectiveness of varying barriers to core vs edge phenotypes
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
	bbMonitorGens=bbMonitorGens)
	
extSeq<-seq(0, 0.5*defBar, 2)
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
			bbMonitorGens)
		save(ext.parameters, extentMat, file="../../Outputs/ExtentTests.RData")
	}
}

#save(ext.parameters, extentMat, file="../../Outputs/ExtentTests.RData")