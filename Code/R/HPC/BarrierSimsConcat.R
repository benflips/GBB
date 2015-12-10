flist<-list.files("../../../../Outputs/", pattern="ExtentTests", full.names=T)

load(flist[1])
nr<-nrow(extentMat)
barSimsMat<-matrix(nrow=length(flist)*nr, ncol=4)

for (ff in 1:length(flist)){
	load(flist[ff])
	rInd<-((ff-1)*nr+1):(ff*nr)
	barSimsMat[rInd,1]<-rep(ext.parameters$defBar, nrow(extentMat))
	barSimsMat[rInd,2]<-rep(ext.parameters$k, nrow(extentMat))
	barSimsMat[rInd,3:4]<-extentMat
}

colnames(barSimsMat)<-c("defBar", "k", "extent", "breachTime")

save(barSimsMat, file="../../../../Outputs/barSimsMat.RData")