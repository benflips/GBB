# walks through parameter combinations and submits jobs to HPC nodes

setwd("/scratch/jc227089/PBSfiles/")
system("rm /scratch/jc227089/PBSfiles/*.*")

Rep<-1:20
dBar<-c(10, 15, 20)
kval<-c(0, 0.1, 0.4)

script.file<-"/home/jc227089/GenBackBurn/GBB/Code/R/HPC/BarrierSimsHPC.R"
fname<-"GBBBarSim"

for (bb in dBar){
	for (ii in Rep) {
		for (kk in kval){
			fid<-paste(bb, "k", kk*10, "rep", ii, sep="")
			##create the sh file
			zz = file(paste(fname, fid,'.sh',sep=''),'w')
			cat('##################################\n',file=zz)
			cat('cd $PBS_O_WORKDIR\n',file=zz)
			cat('source /etc/profile.d/modules.sh\n',file=zz) ### Runs an .sh file which allows modules to be loaded
			cat('module load R\n', file=zz)
			cat("R CMD BATCH --no-save --no-restore '--args Rep=",ii, " defBar=", bb, " k=", kk, "' ", sep="", file=zz)
			cat(script.file, " ", paste(fname, fid,'.Rout',sep=''), "\n", sep="",file=zz)
			cat('##################################\n',file=zz)
			close(zz)
		
			#submit the job
			system(paste('qsub -l nodes=1:ppn=1 -l pmem=1gb -l walltime=36:00:00 ', paste(fname, fid,".sh",sep=""),sep=""))
		}
	}
}