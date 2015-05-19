# Produce images of a GBB over time
rm(list=ls())


source("GBBfunctions.R")
source("GlobalParameters.R")

initGens<-50

run.sim.barr(n, 
			nLoci, 
			eSize, 
			Rmax, 
			Nstar, 
			k, 
			initGens, 
			Ve, 
			mu, 
			defBar, 
			extent=6, 
			lead,
			bbMonitorGens=20,
			save.out="../../Outputs/Demo_Staged.RData")



