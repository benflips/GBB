source("PlotterFunctions.R")
source("GlobalParameters.R")

load("../../Outputs/BaseSpread.RData")

# Plot the trade-off in a final population
plotTradeOff(runList=spreadReps[[1]], 
	parList=parameters, 
	file="../../Figures/TradeOff.pdf")

# Plot of a single realisation of spread.
plotRealisation(runList=spreadReps[[1]], file="../../Figures/Realisation.pdf")

plotRealisation2(runList=spreadReps[[1]], file="../../Figures/Realisation2.pdf")

plotRealisation2_RT(runList=spreadReps[[1]], file="../../Figures/Realisation2_RT.pdf")

plotBasicReps(spreadReps, file="../../Figures/BasicReps.pdf")

source("GBBBaseSpreadStats.R")

plotKernels(poptail, popfront, file="../../Figures/DispersalPhenotypes.pdf")

load("../../Outputs/BarrierTests.RData")

plotVarBarrs(frontMat, coreMat, frontMatEvol, coreMatEvol, file="../../Figures/VarBarrs.pdf", monitorGens)

load("../../Outputs/barSimsMat.RData")

plotBarSims(barSimsMat, file="../../Figures/BarSims.pdf")

load("../../Outputs/Demo_Staged.RData")

plotGBB(pop.out, ints=c(1, 5, 20), file="../../Figures/GBB_staged.pdf")


