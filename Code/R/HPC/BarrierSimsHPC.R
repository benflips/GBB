# HPC Script
# examines the effectiveness of constant barrier to varying introduction extents.
rm(list=ls())
sourc.dir<-""
source(paste(source.dir, "GBBfunctions.R", sep=""))
source(paste(source.dir, "GlobalParameters.R", sep=""))

args=(commandArgs(TRUE))

#evaluate the arguments
# input argument will be Rep
for(i in 1:length(args)) {
	 eval(parse(text=args[[i]]))
}

bbReps<-1

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
		save(ext.parameters, extentMat, 
			file=paste("../../Outputs/ExtentTests", defBar, "_", Rep, ".RData", sep=""))
}
