source("PlotterFunctions.R")

load("../../Outputs/BaseSpread.RData")

# Plot of a single realisation of spread.
plotRealisation(runList=spreadReps[[1]], file="../../Figures/Realisation.pdf")

plotBasicReps(spreadReps, file="../../Figures/BasicReps.pdf")