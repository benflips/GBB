require(fields)
source("GBBfunctions.R")

# Takes a list coming out of run.sim and plots that realisation
plotRealisation<-function(runList, file=NULL){
	if (!is.null(file)) pdf(file=file, height=18, width=9)
	par(mfrow=c(2, 1), mar=c(5,5,2,2), cex.lab=1.4)
	
	X<-runList$pop[,"X"] %/% 3
	Xplot<-tapply(runList$pop[,"X"], X, mean)
	Y<-tapply(dens(runList$pop), X, mean)
	
	plot(Y~Xplot,
		xlab=expression(Distance~from~introduction~(italic(x))),
		ylab="Mean density",
		bty="l", 
		type="l")
	
	
	Y<-tapply(runList$pop[,"di"], X, mean)
	Y2<-tapply(runList$pop[,"di"], X, range)
	Y2<-do.call("rbind", Y2)
	Y2[,2]<-rev(Y2[,2])
	Ypoly<-as.vector(Y2)
	Xpoly<-c(Xplot, rev(Xplot))
	
	plot(Y~Xplot,
		xlab=expression(Distance~from~introduction~(italic(x))),
		ylab=expression(Dispersal~phenotype~(italic(d[i]))),
		bty="l",
		type="l",
		ylim=range(Ypoly))
	polygon(Xpoly, Ypoly,
		col=gray(0.5, 0.5) ,
		border=NA)
		
	if (!is.null(file)) dev.off()
}


# Takes a list coming out of run.sim and plots that realisation
# Modified by RT 05/05/15 to have pretty x axis and color ramp
# Modified further by Ben
plotRealisation2_RT<-function(runList, file=NULL, col.levels=30){
  if (!is.null(file)) pdf(file=file, height=9, width=18)
  par(mar=c(5,5,2,2), cex.lab=1.4)
  
  X<-runList$pop[,"X"] %/% 3
  Xplot<-tapply(runList$pop[,"X"], X, mean)
  Y<-tapply(dens(runList$pop), X, mean)
  Y2<-(tapply(runList$pop[,"di"], X, mean))
  #bpoints<-seq(min(Y2), max(Y2), length.out=31)
  #bpoints[31]<-5
  bpoints<-seq(0, 1, length.out=31)
  bpoints<-quantile(Y2, probs=bpoints)
  Y2cut<-cut((Y2), col.levels, breaks=bpoints, labels=FALSE)
  Y2colRamp<-colorRampPalette(rev(heat.colors(col.levels))) 
  Y2col<-Y2colRamp(col.levels)
  
  bp<-barplot(height=Y,
              xaxt="n",
              xlab=expression(Distance~from~introduction~(italic(x))),
              ylab="Mean density",
              bty="l",
              col=Y2col[as.numeric(Y2cut)],
              border=NA,
              space=0 )
  axis(1,at=round(bp)[ round(bp) %% 10 == 0 ])
  image.plot(smallplot=c(0.25,0.75,0.92,0.94), legend.only=TRUE,zlim=c(0, 100),
             horizontal=T,col= Y2col,legend.line = -2.5,legend.shrink=0.1,
             legend.lab=expression(Dispersal~phenotype~quantile)) 

                
  if (!is.null(file)) dev.off()
}

#takes a timeList and creates a panel at specified intervals
plotGBB<-function(timeList, ints, file=NULL, col.levels=30){
	npan<-length(ints)
	if (!is.null(file)) pdf(file=file, height=9*npan, width=18)
	par(mar=c(5,5,2,2), cex.lab=2, mfrow=c(npan, 1))
	
	X<-timeList[[ints[length(ints)]]][,"X"] %/% 3
	Y2<-(tapply(timeList[[ints[length(ints)]]][,"di"], X, mean))
	bpoints<-seq(min(Y2), max(Y2), length.out=31)
	bpoints[31]<-5
	Y2colRamp<-colorRampPalette(rev(heat.colors(col.levels))) 
	Y2col<-Y2colRamp(col.levels)

  	for (ii in 1:length(ints)){
  		X<-timeList[[ints[ii]]][,"X"] %/% 3
  		Xplot<-tapply(timeList[[ints[ii]]][,"X"], X, mean)
		Y<-tapply(dens(timeList[[ints[ii]]]), X, mean)
		Y2<-(tapply(timeList[[ints[ii]]][,"di"], X, mean))
		Y2cut<-cut((Y2), col.levels, breaks=bpoints, labels=FALSE)
  		
  		xl<-switch((ii==length(ints))+1, "", expression(Distance~from~introduction~(italic(x))))
		bp<-barplot(height=Y,
		  xaxt="n",
		  xlab=xl,
		  ylab="Mean density",
		  bty="l",
		  col=Y2col[as.numeric(Y2cut)],
		  border=NA,
		  space=0,
		  xlim=c(0, (max(Xplot)+30)/3))
		polygon(x=c(max(bp), rep(max(bp)+5, 2), max(bp)),
			y=c(0,0,10,10),
			col=1, 
			border=NA)
		text(bp[30], 50, labels=paste("Generation ", 50+ints[ii]), cex=5)	
		if (ii==length(ints)) {
			bp<-c(bp, max(bp+10))
			axis(1,at=round(bp)[ round(bp) %% 10 == 0 ],
				labels=round(bp)[ round(bp) %% 10 == 0 ]*3)
		}
		if (ii==1) image.plot(smallplot=c(0.25,0.75,0.92,0.94), legend.only=TRUE,zlim=c(min((Y2)),max((Y2))),
			horizontal=T,col= Y2col,legend.line = -2.5,legend.shrink=0.1,
			legend.lab=expression(Mean~dispersal~phenotype~(italic(d[i]))))
	}
	
	if (!is.null(file)) dev.off() 
}


plotTradeOff<-function(runList, parList, file=NULL){
	if (!is.null(file)) pdf(file=file)
	par(mar=c(5,5,2,2), cex.lab=1.4)
	
	X<-seq(0, exp(parList$eSize*parList$nLoci*2), 0.1)
	a<-tOff(parList$Rmax, parList$Nstar, log(X), parList$k)
	densLev<-seq(0, parList$Nstar, parList$Nstar/10)
	Wmat<-c()
	for (dd in 1:length(densLev)){
		Wmat<-cbind(Wmat, bevHolt(densLev[dd], parList$Rmax, parList$Nstar, a))
	}
	matplot(X, Wmat,
		xlab=expression(italic(e^d[i])),
		ylab=expression(E~(italic(W[i]))),
		bty="l", 
		type="l",
		col=gray(seq(0, 1, length=1.3*length(densLev))),
		lty=1,
		lwd=2.5)
		
	if (!is.null(file)) dev.off()
}

# Takes a list of run.sim lists and plots out mean and variation in results
plotBasicReps<-function(repList, file=NULL){
	pop<-c()
	for (ii in 1:length(repList)){
		pop<-rbind(pop, repList[[ii]]$pop)
	}
	
	if (!is.null(file)) pdf(file=file, height=18, width=9)
	par(mfrow=c(2, 1), mar=c(5,5,2,2), cex.lab=1.4)
	
	X<-pop[,"X"] %/% 3
	Xplot<-as.numeric(levels(as.factor(X)))
	Y<-tapply(dens(pop), X, mean)/length(repList)
	#Y2<-tapply(dens(pop), X, range)
	#Y2<-do.call("rbind", Y2)
	#Y2[,2]<-rev(Y2[,2])
	#Ypoly<-as.vector(Y2)
	Xpoly<-c(Xplot, rev(Xplot))
	
	plot(Y~Xplot,
		xlab=expression(Distance~from~introduction~(italic(x))),
		ylab="Mean density",
		bty="l", 
		type="l")
	#polygon(Xpoly, Ypoly,
	#	col=gray(0.5, 0.5) ,
	#	border=NA)
	
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
		
	if (!is.null(file)) dev.off()
}

plotKernels<-function(poptail, popfront, file=NULL){
	dtail<-density(poptail[,"meanD"])
	dfront<-density(popfront[,"meanD"])
	if (!is.null(file)) pdf(file=file, height=9, width=9)
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
	if (!is.null(file)) dev.off()
}

#creates the figure comparing different barriers against core and frontal gene freqs.
plotVarBarrs<-function(frontMat, coreMat, frontMatEvol, coreMatEvol, file=NULL, monitorGens){
	if (!is.null(file)) pdf(file=file, height=18, width=9)
	par(mfrow=c(2,1), mar=c(5,5,2,2), cex.lab=1.4)
	
	Yf<-frontMat[,2]==monitorGens
	Xf<-frontMat[,1]
	fMod<-glm(Yf~Xf, family=binomial)
	Yf.e<-frontMatEvol[,2]==monitorGens
	Xf.e<-frontMatEvol[,1]
	fMod.e<-glm(Yf.e~Xf.e, family=binomial)
	Yc<-coreMat[,2]==monitorGens
	Xc<-coreMat[,1]
	cMod<-glm(Yc~Xc, family=binomial)
	Yc.e<-coreMatEvol[,2]==monitorGens
	Xc.e<-coreMatEvol[,1]
	cMod.e<-glm(Yc.e~Xc.e, family=binomial)
	X<-seq(1,monitorGens,0.2)
	Yf<-monitorGens*predict(fMod, newdata=list("Xf"=X), type="response")
	Yc<-monitorGens*predict(cMod, newdata=list("Xc"=X), type="response")
	Yf.e<-monitorGens*predict(fMod.e, newdata=list("Xf.e"=X), type="response")
	Yc.e<-monitorGens*predict(cMod.e, newdata=list("Xc.e"=X), type="response")
	
	plotMat<-cbind(coreMat, frontMat[,2])
	sunflowerplot(plotMat[,1], plotMat[,2:3],
		pch=1:2,
		col=1:2,
		cex=0.8,
		bty="l",
		xlab="",
		ylab="Barrier strength (generations until breach)")
	lines(X, Yc)
	lines(X, Yf, col=2)
	text(2, 50, labels="A)", cex=1.5)
	legend(0, 45, 
		legend=c("Core", "Front"), 
		title="Gene Frequencies\n(no evolution)", 
		pch=1:2,
		col=1:2, 
		cex=0.8,
		bty="n")
		
	plotMat<-cbind(coreMatEvol, frontMatEvol[,2])
	sunflowerplot(plotMat[,1], plotMat[,2:3],
		pch=1:2,
		col=1:2,
		cex=0.8,
		bty="l",
		xlab="Barrier Size",
		ylab="Barrier strength (generations until breach)")
	lines(X, Yc.e)
	lines(X, Yf.e, col=2)
	text(2, 50, labels="B)", cex=1.5)
	legend(0, 45, 
		legend=c("Core", "Front"), 
		title="Gene Frequencies\n(with evolution)", 
		pch=1:2,
		col=1:2, 
		cex=0.8,
		bty="n")
	
	if (!is.null(file)) dev.off()
}

plotBarSims<-function(bsm, file=NULL){
	if (!is.null(file)) pdf(file=file, height=9, width=9)
	blev<-as.numeric(levels(as.factor(bsm[,"defBar"])))
	klev<-as.numeric(levels(as.factor(bsm[,"k"])))
	colvec<-1:length(blev)
	pchvec<-1:length(blev)
	lwdvec<-1:length(klev)
	par(mar=c(5,5,2,2), cex.lab=1.4)
	plot(bsm[,"extent"], bsm[,"breachTime"]/50,
		type="n",
		bty="l",
		xlab="Backburn extent",
		ylab="Probability of barrier success")
	X<-unique(bsm[,"extent"])
	for (bb in 1:length(blev)){
		cat("bb=",bb,"\n")
		for (kk in 1:length(klev)){
		cat("  kk=",kk, "\n")
			temp<-bsm[(bsm[,"defBar"]==blev[bb] & bsm[,"k"]==klev[kk]),]
			Y<-tapply(temp[,"breachTime"], temp[,"extent"], mean)/50
			points(X, Y, 
				col=bb,
				pch=bb)
			lines(X, Y, 
				col=bb,
				pch=bb,
				lwd=kk)
		}
	}
	legend('right', 
		pch=pchvec,
		col=colvec,
		legend=blev,
		title="Barrier width",
		bty="n")
	if (!is.null(file)) dev.off()
}