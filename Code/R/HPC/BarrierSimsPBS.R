setwd("/scratch/jc227089/PBSfiles/")


R0<-seq(6, 10, 2)
Nstar<-c(10)
pop.reps<-20

script.file<-"/scratch/jc227089/evo-dispersal/evostoch/2D_alphawalker.R"
fname<-"evostoch2D"
for (kk in 1: length(Nstar)){
	for (jj in 1:length(R0)){
		for (ii in 1:pop.reps) {

			fid<-paste(R0[jj], Nstar[kk], ii, sep="_")
			##create the sh file
			zz = file(paste(fname, fid,'.sh',sep=''),'w')
			cat('##################################\n',file=zz)
			cat('cd $PBS_O_WORKDIR\n',file=zz)
			cat('source /etc/profile.d/modules.sh\n',file=zz) ### Runs an .sh file which allows modules to be loaded
			cat('module load R\n', file=zz)
			cat("R CMD BATCH --no-save --no-restore '--args repn=",ii, " R0=",R0[jj], " Nstar=",Nstar[kk], "' ", sep="", file=zz)
			cat(script.file, " ", paste(fname, fid,'.Rout',sep=''), "\n", sep="",file=zz)
			cat('##################################\n',file=zz)
			close(zz)
				
			#submit the job
			system(paste("qsub -m n ", paste(fname, fid,".sh",sep=""),sep=""))
		}
	}
}
