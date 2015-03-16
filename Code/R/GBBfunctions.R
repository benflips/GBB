# Initialises the simulation
init<-function(n, nLoci){
	X<-rep(0, n)
	alleles<-matrix(rbinom(n*nLoci, 1, 0.5), nrow=n,
		dimnames=list(NULL, paste("a", 1:nLoci, sep="")))
	cbind(X, alleles)
}

# Calculates dispersal phenotype based on individual alleles, their effect, and Ve
dPhen<-function(pop, eSize, Ve){
	gtypes<-pop[,grepl("a", colnames(pop))]
	gtypes<-apply(gtypes, 1, sum)*eSize
	ptypes<-rnorm(length(gtypes), mean=gtypes, sd=sqrt(Ve))
	ptypes
}

# calculates fitness given Rmax and phenotype
tOff<-function(Rmax, phen, k=1){
	Rmax*exp(-k*phen)
}

run.sim<-function(n, nLoci, eSize, Ve, Rmax){
	pop<-init(n, nLoci)
	disp<-dPhen(pop, eSize, Ve)
	fitness<-tOff(Rmax, disp)
	list(x=disp, y=fitness)
}

