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
tOff<-function(phen, k=1){
	exp(-k*phen)
}

# Beverton Holt population growth
bevHolt<-function(N, Rmax, Nstar){
  a<-(Rmax-1)/Nstar
  Rmax/(1+a*N)
}

# reproduces individuals
reproduce<-function(pop){
	
}



run.sim<-function(n, nLoci, eSize, Rmax){
	init(n, nLoci, eSize)
	#disp<-dPhen(pop, eSize, Ve)
	#fitness<-tOff(disp)
	#Eo<-bevHolt(nrow(pop), )
	#list(x=disp, y=fitness)
}

