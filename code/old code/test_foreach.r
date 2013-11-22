
rm(list=ls())

	
library(foreach)
library(doSNOW)


dir.create("results2")

workerlist <- rep("localhost", times = 23)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)


foreach(iteration = 1:46) %dopar% {
#
temp <- rnorm(20)
save(temp, file = paste("results2/temp", iteration, sep = "_"))
	

gc() 
}



stopCluster(cl)

