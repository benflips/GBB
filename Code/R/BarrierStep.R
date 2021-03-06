# examines the effectiveness of varying barriers to core vs edge phenotypes
rm(list=ls())
source("GBBfunctions.R")
source("GlobalParameters.R")
load(file="../../Outputs/BaseSpreadSummary.RData")

lead<-5
maxX<-30

parameters<-list(n=n, 
	nLoci=nLoci, 
	eSize=eSize,
	Rmax=Rmax,
	Nstar=Nstar,
	k=k,
	Ve=Ve,
	mu=mu,
	maxX=maxX,
	lead=lead,
	initGens=initGens,
	monitorGens=monitorGens)
	


barSeq<-minBarSize:maxBarSize

# Test frontal genefreqs against varying barriers

frontFreqs<-popfront[,2:(2*nLoci+1)]
frontFreqs<-apply(frontFreqs, 2, function(x){sum(x)/length(x)})

frontMat<-cbind(barSize=rep(barSeq, times=rep(baseReps, length(barSeq))), breachTime=rep(NA, length(barSeq)*baseReps))
cat("Frontal genotypes, no evolution...", "\n")
for (bb in 1:length(barSeq)){
	cat("Working through barrier", bb, "of", length(barSeq), "...\n")
	for (rr in 1:baseReps){
	cat("Replicate", rr, "\n")
		frontMat[((bb-1)*baseReps+rr),2]<-bar.test.sim(n*maxX, 
			nLoci, 
			eSize, 
			Rmax, 
			Nstar, 
			Ve, 
			mu, 
			gFreqs=frontFreqs, 
			maxX=maxX, 
			barSize=barSeq[bb],
			lead=lead, 
			monitorGens)
	}
}

frontMatEvol<-cbind(barSize=rep(barSeq, times=rep(baseReps, length(barSeq))), breachTime=rep(NA, length(barSeq)*baseReps))
cat("Frontal genotypes, with evolution...", "\n")
for (bb in 1:length(barSeq)){
	cat("Working through barrier", bb, "of", length(barSeq), "...\n")
	for (rr in 1:baseReps){
	cat("Replicate", rr, "\n")
		frontMatEvol[((bb-1)*baseReps+rr),2]<-bar.test.sim.evol(
			n*maxX, 
			nLoci, 
			eSize, 
			Rmax, 
			Nstar,
			k, 
			Ve, 
			mu, 
			gFreqs=frontFreqs, 
			maxX=maxX, 
			barSize=barSeq[bb],
			lead=lead, 
			monitorGens)
	}
}



# Test core genefreqs against varying barriers

coreFreqs<-poptail[,2:(2*nLoci+1)]
coreFreqs<-apply(coreFreqs, 2, function(x){sum(x)/length(x)})

coreMat<-cbind(barSize=rep(barSeq, times=rep(baseReps, length(barSeq))), breachTime=rep(NA, length(barSeq)*baseReps))
cat("Core genotypes, no evolution...", "\n")
for (bb in 1:length(barSeq)){
	cat("Working through barrier", bb, "of", length(barSeq), "...\n")
	for (rr in 1:baseReps){
	cat("Replicate", rr, "\n")
		coreMat[((bb-1)*baseReps+rr),2]<-bar.test.sim(n*maxX, 
			nLoci, 
			eSize, 
			Rmax, 
			Nstar, 
			Ve, 
			mu, 
			gFreqs=coreFreqs, 
			maxX=maxX, 
			barSize=barSeq[bb],
			lead=lead, 
			monitorGens)
	}
}

coreMatEvol<-cbind(barSize=rep(barSeq, times=rep(baseReps, length(barSeq))), breachTime=rep(NA, length(barSeq)*baseReps))
cat("Core genotypes, with evolution...", "\n")
for (bb in 1:length(barSeq)){
	cat("Working through barrier", bb, "of", length(barSeq), "...\n")
	for (rr in 1:baseReps){
	cat("Replicate", rr, "\n")
		coreMatEvol[((bb-1)*baseReps+rr),2]<-bar.test.sim.evol(
			n*maxX, 
			nLoci, 
			eSize, 
			Rmax, 
			Nstar,
			k, 
			Ve, 
			mu, 
			gFreqs=coreFreqs, 
			maxX=maxX, 
			barSize=barSeq[bb],
			lead=lead, 
			monitorGens)
	}
}

save(parameters, frontMat, frontMatEvol, coreMat, coreMatEvol, file="../../Outputs/BarrierTests.RData")