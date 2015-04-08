source("PlotterFunctions.R")

load("../../Outputs/BaseSpread.RData")

# Plot the trade-off in a final population
plotTradeOff(runList=spreadReps[[1]], 
	parList=parameters, 
	file="../../Figures/TradeOff.pdf")

# Plot of a single realisation of spread.
plotRealisation(runList=spreadReps[[1]], file="../../Figures/Realisation.pdf")

plotBasicReps(spreadReps, file="../../Figures/BasicReps.pdf")