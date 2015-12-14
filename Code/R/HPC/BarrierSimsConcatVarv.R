flist<-list.files("../../../../Outputs/", pattern="ExtentTests", full.names=T)

vVec<-substr(flist, regexpr("v", flist)+1, regexpr("v", flist)+1)
vVec<-as.numeric(ifelse(vVec=="1", "10", vVec))
load(flist[1])
nr<-nrow(extentMat)
barSimsMat<-matrix(nrow=length(flist)*nr, ncol=4)

for (ff in 1:length(flist)){
	load(flist[ff])
	rInd<-((ff-1)*nr+1):(ff*nr)
	barSimsMat[rInd,1]<-rep(ext.parameters$defBar, nrow(extentMat))
	barSimsMat[rInd,2]<-rep(vVec[ff], nrow(extentMat))
	barSimsMat[rInd,3:4]<-extentMat
}

colnames(barSimsMat)<-c("defBar", "v", "extent", "breachTime")

barSimsMatVarv<-barSimsMat
load(file="../../../../Outputs/barSimsMat.RData")
barSimsMat<-barSimsMat[barSimsMat[,"k"]==0.2, ]
barSimsMat[,"k"]<-Inf
barSimsMatVarv<-rbind(barSimsMatVarv, barSimsMat)

save(barSimsMatVarv, file="../../../../Outputs/barSimsMatVarv.RData")