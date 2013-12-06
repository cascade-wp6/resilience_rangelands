
rm(list=ls())

	
library(foreach)
library(doSNOW)


dir.create("results_test")

winWorker <- list(host = "localhost", 
				rscript = "C:\\R\\R-2.15.3\\bin\\x64\\Rscript.exe", 
				snowlib = "C:\\R\\R-2.15.3\\library"
	)
winWorker <- list(host = "localhost", 
				rscript = "C:/R/R-2.15.3/bin/x64/Rscript.exe", 
				snowlib = "C:/R/R-2.15.3/library"
	)

workstation <-  list(host = "162.38.184.118", user = "schneider",
         rscript = "/usr/lib/R/bin/Rscript",
		 snowlib = "/usr/lib/R/library/")
		
workerlist <- rep(list(workstation), times = 20)
workerlist <- c( rep(list(winWorker), times = 10),rep(list(workstation), times = 23))

#workerlist <- rep("localhost", times = 23)

cl <- makeSOCKcluster(workerlist, master = "162.38.184.88")

clusterCall(cl,function() Sys.info()[c("nodename","machine")])

# /usr/lib/R/bin/Rscript /usr/lib/R/library//snow/RSOCKnode.R  MASTER=LAB PORT=10187 OUT=/dev/null SNOWLIB=/usr/lib/R/library/ 
registerDoSNOW(cl)


foreach(iteration = 1:46, .combine = rbind) %dopar% {
#
temp <- rnorm(10)

return(temp)
gc() 
}



stopCluster(cl)

