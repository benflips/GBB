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
tOff<-function(phen, k=0.1){
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
reproduce<-function(pop, Rmax, Nstar, nLoci){
#browser()
	EW<-bevHolt(nrow(pop), Rmax, Nstar)*tOff(pop[,"di"])
	W<-rpois(length(EW), EW)
	gtype.loc<-grepl("a", colnames(pop))
	pop2<-c()
	for (ii in 1:nrow(pop)){
		for (ww in 1:W[ii]){
			temp<-pop[ii,]
			tempF<-temp
			tempM<-pop[sample(1:nrow(pop), 1),]
			temp[gtype.loc]<-as.vector(rbind(gametes(tempF), gametes(tempM)))
			pop2<-rbind(pop2, temp)
		}
	}
	pop2
}


run.sim<-function(n, nLoci, eSize, Rmax, Nstar){
	temp<-init(n, nLoci, eSize)
	pop<-temp$pop
	Ve<-temp$Ve
	rm(temp)
	pop2<-reproduce(pop, Rmax, Nstar, nLoci)
	list(pop, pop2)
}

