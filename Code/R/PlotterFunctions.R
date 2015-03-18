source("GBBfunctions.R")

# Takes a list coming out of run.sim and plots that realisation
plotRealisation<-function(runList, file){
	pdf(file=file, height=18, width=9)
	par(mfrow=c(2, 1), mar=c(5,5,2,2), cex.lab=1.4)
	
	X<-runList$pop[,"X"] %/% 3
	Xplot<-as.numeric(levels(as.factor(X)))
	Y<-tapply(dens(runList$pop), X, mean)
	
	plot(Y~Xplot,
		xlab=expression(Distance~from~introduction~(italic(x))),
		ylab="Mean density",
		bty="l", 
		type="l")
	
	Y<-tapply(runList$pop[,"di"], X, mean)
	plot(Y~Xplot,
		xlab=expression(Distance~from~introduction~(italic(x))),
		ylab=expression(Mean~dispersal~phenotype~(italic(d[i]))),
		bty="l",
		type="l")
		
	dev.off()
}

# Takes a list of run.sim lists and plots out mean and variation in results
plotBasicReps<-function(repList, file){
	pop<-c()
	for (ii in 1:length(repList)){
		pop<-rbind(pop, repList[[ii]]$pop)
	}
	
	pdf(file=file, height=18, width=9)
	par(mfrow=c(2, 1), mar=c(5,5,2,2), cex.lab=1.4)
	
	X<-pop[,"X"] %/% 3
	Xplot<-as.numeric(levels(as.factor(X)))
	Y<-tapply(dens(pop), X, mean)
	Y2<-tapply(dens(pop), X, range)
	Y2<-do.call("rbind", Y2)
	Y2[,2]<-rev(Y2[,2])
	Ypoly<-as.vector(Y2)
	Xpoly<-c(Xplot, rev(Xplot))
	
	plot(Y~Xplot,
		xlab=expression(Distance~from~introduction~(italic(x))),
		ylab="Mean density",
		bty="l", 
		type="l", 
		ylim=range(Ypoly))
	polygon(Xpoly, Ypoly,
		col=gray(0.5, 0.5) ,
		border=NA)
	
	Y<-tapply(pop[,"di"], X, mean)
	Y2<-tapply(pop[,"di"], X, range)
	Y2<-do.call("rbind", Y2)
	Y2[,2]<-rev(Y2[,2])
	Ypoly<-as.vector(Y2)
	
	plot(Y~Xplot,
		xlab=expression(Distance~from~introduction~(italic(x))),
		ylab=expression(Mean~dispersal~phenotype~(italic(d[i]))),
		bty="l",
		type="l",
		ylim=range(Ypoly))
	polygon(Xpoly, Ypoly,
		col=gray(0.5, 0.5) ,
		border=NA)
		
	dev.off()
}

