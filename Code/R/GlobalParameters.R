# define global parameters

## Demographic and evolutionary pars
nLoci<-20
initFreq<-0.5
Rmax<-5
Nstar<-200
n<-Nstar
mu<-log(2) # log mean dispersal distance (at init) 
VT<-log(1.5)^2 # Total variance in di
h2<-0.99
Vg<-h2*VT
Ve<-(1-h2)*VT
eSize<-sqrt((Vg)/((2*nLoci)*initFreq*(1-initFreq)))
k<-0.2


## Experiment descriptors
initGens<-50
monitorGens<-50
baseReps<-1
bbReps<-20
maxBarSize<-10

