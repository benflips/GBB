#source("Plots.R")

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

dtail<-density(poptail[,"meanD"])
dfront<-density(popfront[,"meanD"])


#pdf(file="../../Figures/DispersalPhenotypes.pdf")
	par(mar=c(5,5,2,2), cex.lab=1.4)
	plot(c(dtail$x, dfront$x), c(dtail$y, dfront$y), 
		type="n",
		xlab="Mean dispersal distance",
		ylab="Probability density",
		bty="l")
	lines(dtail, lty=1)
	lines(dfront, lty=2)
	legend('topright', 
		legend=c("Range core", "Expanding front"), 
		lty=1:2,
		bty="n")
#dev.off()