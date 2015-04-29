# define global parameters

## Demographic and evolutionary pars
nLoci<-20
initFreq<-0.5
Rmax<-5
Nstar<-200
n<-Nstar
mu<-log(2) # log mean dispersal distance (at init) 
VT<-log(1.5)^2 # Total variance in di
h2<-0.3
Vg<-h2*VT
Ve<-(1-h2)*VT
eSize<-sqrt((Vg)/((2*nLoci)*initFreq*(1-initFreq)))
k<-0.2


## Experiment descriptors
initGens<-10 #initial spread generations
monitorGens<-50 # monitoring generations for BarrierStep
baseReps<-2 # basic number of replicates of each simulation
bbReps<-20 # replicates of backburn simulations
maxBarSize<-50
minBarSize<-1
defBar<-20 # barrier size for constant barrier sims
lead<-10 # how far in front of population you set up the barrier/burn
bbMonitorGens<-50 # monitor generations for backburn sims

