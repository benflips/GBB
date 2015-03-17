# Initialises the simulation
init<-function(n, nLoci, eSize){
	X<-runif(n, 0, 1)
	alleles<-matrix(rbinom(2*n*nLoci, 1, 0.5), nrow=n,
		dimnames=list(NULL, paste("a", 1:(2*nLoci), sep="")))
	pop<-cbind(X, alleles, di=rep(0, n))
	# work out appropriate Ve such that h2 ca 0.3
	gtypes<-pop[,grepl("a", colnames(pop))]
	gtypes<-apply(gtypes, 1, sum)*eSize
	Vg<-var(gtypes)
	Ve<-Vg/0.3
	pop<-dPhen(pop, eSize, Ve)
	list(pop=pop, Ve=Ve)
}

# Calculates dispersal phenotype based on individual alleles, their effect, and Ve
dPhen<-function(pop, eSize, Ve){
	gtypes<-pop[,grepl("a", colnames(pop))]
	gtypes<-apply(gtypes, 1, sum)*eSize
	ptypes<-rnorm(length(gtypes), mean=gtypes, sd=sqrt(Ve))
	pop[,"di"]<-exp(ptypes)
	pop
}

# calculates relative fitness (max=1) given dispersal phenotype
tOff<-function(phen, k=0.01){
	exp(-k*phen)
}

# Beverton Holt population growth
bevHolt<-function(N, Rmax, Nstar){
  a<-(Rmax-1)/Nstar
  Rmax/(1+a*N)
}

# Generates gametes for an individual
gametes<-function(poprow){
	gtypes<-poprow[grepl("a", names(poprow))]
	temp<-rbinom(nLoci, 1, 0.5)
	temp<-as.vector(rbind(temp, temp!=1))
	gametes<-gtypes[as.logical(temp)]
	gametes
}

# reproduces individuals
reproduce<-function(pop, Rmax, Nstar, nLoci, eSize, Ve){
	EW<-bevHolt(dens(pop), Rmax, Nstar)*tOff(pop[,"di"])
	W<-rpois(length(EW), EW)
	gtype.loc<-grepl("a", colnames(pop))
	pop2<-c()
	for (ii in 1:nrow(pop)){
		if (W[ii]==0) next
		for (ww in 1:W[ii]){
			temp<-pop[ii,]
			tempF<-temp
			popSet<-which((pop[,"X"] %/% 1) == (tempF["X"] %/% 1))
			tempM<-pop[sample(popSet, 1),]
			temp[gtype.loc]<-as.vector(rbind(gametes(tempF), gametes(tempM)))
			pop2<-rbind(pop2, temp)
		}
	}
	pop<-dPhen(pop2, eSize, Ve)
}

# moves individuals on the lattice
disperse<-function(pop){
	pop[,"X"]<-rnorm(nrow(pop), mean=pop[,"X"], sd=pop[,"di"])
	pop[,"X"]<-ifelse(pop[,"X"]<0, -pop[,"X"], pop[,"X"])
	pop
}

# calculates density back to the individual
dens<-function(pop){
	Xbin<-(pop[,"X"] %/% 1)+1
	dens<-table(factor(Xbin, levels=1:max(Xbin)))
	dens[Xbin]
}

run.sim<-function(n, nLoci, eSize, Rmax, Nstar){
	temp<-init(n, nLoci, eSize)
	pop<-temp$pop
	Ve<-temp$Ve
	rm(temp)
	N<-nrow(pop)
	mean.di<-mean(pop[,"di"])
	X.g<-var(pop[,"X"])
	X.min<-min(pop[,"X"])
	for (gg in 1:100){
		cat("Working through generation ", gg, "\n")
		pop<-reproduce(pop, Rmax, Nstar, nLoci, eSize, Ve)
		pop<-disperse(pop)
		N<-c(N, nrow(pop))
		mean.di<-c(mean.di, mean(pop[,"di"]))
		X.g<-c(X.g, var(pop[,"X"]))
		X.min<-c(X.min, min(pop[,"X"]))
	}
	list(N=N, mean.di=mean.di, X.g=X.g, X.min=X.min, pop=pop)
}

