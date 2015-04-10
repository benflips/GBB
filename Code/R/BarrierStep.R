# examines the effectiveness of varying barriers to core vs edge phenotypes
rm(list=ls())
source("GBBfunctions.R")
source("GlobalParameters.R")

parameters<-list(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	initGens=initGens)
	
load(file="../../Outputs/BaseSpreadSummary.RData")

barSeq<-minBarSize:maxBarSize

# Test frontal genefreqs against varying barriers

frontFreqs<-popfront[,2:(2*nLoci+1)]
frontFreqs<-apply(frontFreqs, 2, function(x){sum(x)/length(x)})

frontMat<-cbind(barSize=rep(barSeq, times=rep(baseReps, length(barSeq))), breachTime=rep(NA, length(barSeq)*baseReps))
for (bb in 1:length(barSeq)){
	cat("Working through barrier", bb, "of", length(barSeq), "...\n")
	for (rr in 1:baseReps){
	cat("Replicate", rr, "\n")
		frontMat[((bb-1)*baseReps+rr),2]<-bar.test.sim(n*10, 
			nLoci, 
			eSize, 
			Rmax, 
			Nstar, 
			Ve, 
			mu, 
			gFreqs=frontFreqs, 
			maxX=10, 
			barSize=barSeq[bb], 
			monitorGens)
	}
}

# Test core genefreqs against varying barriers

coreFreqs<-poptail[,2:(2*nLoci+1)]
coreFreqs<-apply(coreFreqs, 2, function(x){sum(x)/length(x)})

coreMat<-cbind(barSize=rep(barSeq, times=rep(baseReps, length(barSeq))), breachTime=rep(NA, length(barSeq)*baseReps))
for (bb in 1:length(barSeq)){
	cat("Working through barrier", bb, "of", length(barSeq), "...\n")
	for (rr in 1:baseReps){
	cat("Replicate", rr, "\n")
		coreMat[((bb-1)*baseReps+rr),2]<-bar.test.sim(n*10, 
			nLoci, 
			eSize, 
			Rmax, 
			Nstar, 
			Ve, 
			mu, 
			gFreqs=coreFreqs, 
			maxX=10, 
			barSize=barSeq[bb], 
			monitorGens)
	}
}

save(frontMat, coreMat, file="../../Outputs/BarrierTests.RData")