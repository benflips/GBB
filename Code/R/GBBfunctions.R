# Initialises the simulation
init<-function(n, nLoci, eSize, Ve, mu){
	X<-runif(n, 0, 1)
	alleles<-matrix(rbinom(2*n*nLoci, 1, 0.5), nrow=n,
		dimnames=list(NULL, paste("a", 1:(2*nLoci), sep="")))
	pop<-cbind(X, alleles, di=rep(0, n))
	pop<-dPhen(pop, eSize, Ve, nLoci, mu)
	pop
}

# Calculates dispersal phenotype based on individual alleles, their effect, and Ve
dPhen<-function(pop, eSize, Ve, nLoci, mu){
	gtypes<-pop[,grepl("a", colnames(pop))]
	gtypes<-apply(gtypes, 1, sum)*eSize-(nLoci*2*0.5*eSize)+mu
	ptypes<-rnorm(length(gtypes), mean=gtypes, sd=sqrt(Ve))
	pop[,"di"]<-ptypes
	pop
}

# returns alpha of the BevHolt function
tOff<-function(Rmax, Nstar, di, k){
	(Rmax-1)/(Nstar*exp(-k*exp(di)))
}

# Beverton Holt population growth
bevHolt<-function(N, Rmax, Nstar, a){
  #a<-(Rmax-1)/Nstar
  Rmax/(1+a*N)
}

# Generates gametes for an individual
gametes<-function(poprow){
	orderlist<-c(seq(1, 2*nLoci, 2), seq(2, 2*nLoci, 2))
	gtypes<-poprow[grepl("a", names(poprow))]
	temp<-rbinom(nLoci, 1, 0.5)
	temp<-c(temp, temp!=1)[order(orderlist)]
	gametes<-gtypes[temp==1]
	gametes
}


# deprecated reproduction function. Newer version is substantially faster at large N
reproduce<-function(pop, Rmax, Nstar, nLoci, eSize, Ve, k, mu){
        orderlist<-order(c(seq(1, 2*nLoci, 2), seq(2, 2*nLoci, 2)))
        alpha<-tOff(Rmax, Nstar,pop[,"di"], k)
        EW<-bevHolt(dens(pop), Rmax, Nstar, alpha)
        W<-rpois(length(EW), EW)
        gtype.loc<-grepl("a", colnames(pop))
       	pop2<-matrix(0, nrow=sum(W), ncol=ncol(pop), dimnames=dimnames(pop))
       	offInd<-1
        for (ii in 1:nrow(pop)){
                if (W[ii]==0) next
                for (ww in 1:W[ii]){
                        temp<-pop[ii,]
                        tempF<-temp
                      	popSet<-which((pop[,"X"] %/% 1) == (tempF["X"] %/% 1))
                        tempM<-pop[sample(popSet, 1),]
                        temp[gtype.loc]<-c(gametes(tempF), gametes(tempM))[orderlist]
                       	pop2[offInd,]<-temp
                       	offInd<-offInd+1
                }
        }
        pop<-dPhen(pop2, eSize, Ve, nLoci, mu)
        pop
 }


# moves individuals on the lattice
disperse<-function(pop){
	pop[,"X"]<-rnorm(nrow(pop), mean=pop[,"X"], sd=exp(pop[,"di"]))
	pop[,"X"]<-ifelse(pop[,"X"]<0, -pop[,"X"], pop[,"X"])
	pop
}

# calculates density back to the individual
dens<-function(pop){
	Xbin<-(pop[,"X"] %/% 1)+1
	dens<-table(factor(Xbin, levels=1:max(Xbin)))
	dens[Xbin]
}

run.sim<-function(n, nLoci, eSize, Rmax, Nstar, k, initGens, Ve, mu){
	pop<-init(n, nLoci, eSize, Ve, mu)
	N<-nrow(pop)
	mean.di<-mean(pop[,"di"])
	X.g<-var(pop[,"X"])
	X.max<-max(pop[,"X"])
	for (gg in 1:initGens){
		cat("Working through generation ", gg, "\n")
		pop<-reproduce(pop, Rmax, Nstar, nLoci, eSize, Ve, k, mu)
		pop<-disperse(pop)
		N<-c(N, nrow(pop))
		mean.di<-c(mean.di, mean(pop[,"di"]))
		X.g<-c(X.g, var(pop[,"X"]))
		X.max<-c(X.max, max(pop[,"X"]))
	}
	list(N=N, mean.di=mean.di, X.g=X.g, X.max=X.max, pop=pop)
}

