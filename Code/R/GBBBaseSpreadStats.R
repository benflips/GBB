# get trailing and leading 1000 individuals from each replicate population
poptail<-c()
popfront<-c()

for (ii in 1:length(spreadReps)){
	pop<-spreadReps[[ii]]$pop
	if (is.null(pop)) next
	pop<-pop[order(pop[,"X"]),]
	poptail<-rbind(poptail, head(pop, 1000))
	popfront<-rbind(popfront, tail(pop, 1000))
}

poptail<-cbind(poptail, meanD=exp(poptail[,"di"]))
popfront<-cbind(popfront, meanD=exp(popfront[,"di"]))
save(poptail, popfront, file="../../Outputs/BaseSpreadSummary.RData")