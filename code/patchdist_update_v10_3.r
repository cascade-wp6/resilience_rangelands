

rm(list=ls())

#setwd("E:\\Eigene Dokumente\\Uni\\projects\\2013 CASCADE grazing\\sim9\\")
setwd("C:\\Users\\SCHNEIDER\\Documents\\projects\\CAS01_grazing\\data\\sim10\\")



# time and resolution of simulation
timesteps = 2000
delta = 1/5

snapshots <-  (timesteps-20:0*50)/delta+1

output <- read.csv("output.csv")


filenames <- list.files("results/")



library(foreach)
library(doSNOW)

#workstation <-  list(host = "162.38.184.118", user = "schneider",
         rscript = "/usr/lib/R/bin/Rscript",
		 snowlib = "/usr/lib/R/library/")
		
#workerlist <- rep(list(workstation), times = 23)

workerlist <-  rep("localhost", times = 10)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)


foreach(i = filenames) %dopar% {
#for(i in filenames) {
#i =   2603  #which(! 1:2640  %in% read.csv("output_fit.csv")$ID)

# output$ID[which(output$rho_plus > 0.4  & output$rho_plus < 0.5)] 
#i = 1405
# output$ID[which(output$rho_plus > 0.68  & output$rho_plus < 0.7)] 
#i =  25508
# output$ID[which(output$rho_plus > 0.02  & output$rho_plus < 0.1)] 
#i = 900# 328 # 2220 # 
#i = 29880 #5053 #
#i = 9266
load(paste("results\\", i, sep = ""))
#load(paste("results\\result_", 17096, sep = ""))



result$fit <- list()
#result$fit$PL <- data.frame(snapshot = NA, a = NA, alpha = NA, third = NA, p_a = NA, p_alpha = NA, p_third = NA, AIC = NA)[-1,]
#result$fit$TPLup <- data.frame(snapshot = NA, a = NA, alpha = NA, b = NA, p_a = NA, p_alpha = NA, p_b = NA, AIC = NA)[-1,]
#result$fit$TPLdown <- data.frame(snapshot = NA, a = NA, alpha = NA, Sx = NA, p_a = NA, p_alpha = NA, p_Sx = NA, AIC = NA)[-1,]
#result$fit$EXP <- data.frame(snapshot = NA, a = NA, eps = NA, third = NA, p_a = NA, p_alpha = NA, p_third = NA,  AIC = NA)[-1,]
#result$fit$best <- vector()

#	for(j in 2:22 ) {
	
	#j = 6
		dd4 <- do.call("rbind", result$cumpatch)

				flag_full = FALSE

				if(result$out$rho_plus < 0.02 ) {
					flag_desert = TRUE
					if(result$out$rho_plus > 0 & !is.na(result$cumpatch[[22]])) {
						lpatches = sapply(2:22, function(x) max(result$cumpatch[[x]]$size))
					} else {
						lpatches = 0 
					}
				} else { 	# excluding deserts from model fit
					dd4 <- do.call("rbind", result$cumpatch)
					b = mean(sapply(2:22, function(x) 1/sum(result$cumpatch[[x]]$n) ))
					lpatches =  sapply(2:22, function(x) max(result$cumpatch[[x]]$size))
					flag_desert = FALSE
					if(mean(sapply(2:22, function(x) length(result$cumpatch[[x]]$size))) < 3) {flag_full = TRUE} else {flag_full = FALSE}
				}
				
				flag = !flag_desert &  !flag_full
						
if(flag) {


	#if(!is.na(result$cumpatch[[j]][1,1])) dd3 <- data.frame(size = result$cumpatch[[j]]$size, n = result$cumpatch[[j]]$p)


	result$fit$AIC <- vector()
	result$fit$summary <- list()
	#############
PLlm <- lm(I(log(p)) ~  I(log(size)) , data = dd4) 

#if( max(dd3$size) <= 7000) {
try({result$fit$PL <- nls(I(log(p)) ~ log(a) - alpha * log(size), 
		data = dd4,
		start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]),
		trace = FALSE,
		#algorithm = "port",
		nls.control(maxiter = 50)
		)}, silent = TRUE
	)
#	} 
	
	if(!is.null(result$fit$PL)) {
		result$fit$AIC[1] <- AIC(result$fit$PL)
		result$fit$summary$PL <- data.frame(a = coefficients(result$fit$PL)[1], alpha = coefficients(result$fit$PL)[2], third = NA, p_a = summary(result$fit$PL)$coefficients[1,4], p_alpha = summary(result$fit$PL)$coefficients[2,4], p_third = NA, AIC = AIC(result$fit$PL))
	} else {
		result$fit$summary$PL  <- list(NA)
		result$fit$PL  <- list(NA)
		result$fit$AIC[1] <- NA
	}

###########

#b=result$cumpatch[[j]]$p[dim(result$cumpatch[[j]])[1]] #1/sum(result$cumpatch[[j]]$n)  
try({result$fit$TPLup <- nls(I(log(p)) ~ I( log(a) - alpha * log(size) + log(1+b/(a*size^(-alpha))) ), 
		data = dd4,
		start = list(a =  exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]) , #, b = 1/sum(result$cumpatch[[j]]$n)
        trace = FALSE,
		#algorithm = "port",
		#lower = c(0, 0), upper = c(1, NA),
		nls.control(maxiter = 50)
		)}, silent = TRUE
	)
	
	
	if(!is.null(result$fit$TPLup)) {
		result$fit$AIC[2] <- AIC(result$fit$TPLup) 
		result$fit$summary$TPLup <- data.frame(a = coefficients(result$fit$TPLup)[1], alpha = coefficients(result$fit$TPLup)[2], b = b, p_a = summary(result$fit$TPLup)$coefficients[1,4], p_alpha = summary(result$fit$TPLup)$coefficients[2,4], p_b = NA, AIC = AIC(result$fit$TPLup))

	} else { 
		result$fit$summary$TPLup  <- list(NA)
		result$fit$TPLup  <- list(NA)
		result$fit$AIC[2] <- NA
	}
	
	
try( {result$fit$TPLdown <- nls(I(log(p)) ~ I( log(a) - alpha * log(size) - (size * Sx) ), 
		data = dd4,
		start = list(a = exp(PLlm$coefficients[1]), alpha =  exp(PLlm$coefficients[2]), Sx = 1/1000),
        #algorithm = "port",
		trace = FALSE
		)}, silent = TRUE
	)		

	if(!is.null(result$fit$TPLdown) & !coefficients(result$fit$TPLdown)[[3]] <= 0) {
		result$fit$AIC[3] <- AIC(result$fit$TPLdown) 
		result$fit$summary$TPLdown <- data.frame(a = coefficients(result$fit$TPLdown)[1], alpha = coefficients(result$fit$TPLdown)[2], Sx = coefficients(result$fit$TPLdown)[3], p_a = summary(result$fit$TPLdown)$coefficients[1,4], p_alpha = summary(result$fit$TPLdown)$coefficients[2,4], p_Sx = summary(result$fit$TPLdown)$coefficients[3,4], AIC = AIC(result$fit$TPLdown))

	} else {
		result$fit$summary$TPLdown <- list(NA)	
		result$fit$TPLdown <- list(NA)
		result$fit$AIC[3] <- NA
	}

###########
	
try( {result$fit$EXP <- nls(I(log(p)) ~ I(log(a) -(eps*size)) , 
		data = dd4,
		start = list(a = exp(PLlm$coefficients[1]) ,eps = 1),
        #algorithm = "port",
		trace = FALSE
		)}, silent = TRUE
	)
	
		
	if(!is.null(result$fit$EXP)) {
		result$fit$AIC[4] <- AIC(result$fit$EXP) 
		result$fit$summary$EXP <- data.frame(a = coefficients(result$fit$EXP)[1], eps = coefficients(result$fit$EXP)[2], third = NA, p_a = summary(result$fit$EXP)$coefficients[1,4], p_alpha = summary(result$fit$EXP)$coefficients[2,4], p_Sx = NA, AIC = AIC(result$fit$EXP))

	} else {
		result$fit$summary$EXP <- list(NA)
		result$fit$EXP <- list(NA)
		result$fit$AIC[4] <- NA
	}
	
		
	#min(models$AIC, na.rm = TRUE)
	result$fit$dAIC <- 	result$fit$AIC -min(result$fit$AIC, na.rm = TRUE)
	
	
	result$fit$best <- which.min(result$fit$AIC[-4]+c(+0,0,0))

	} else {
	
		if(flag_desert) { result$fit$best <- 0 } 
		if(flag_full)	{result$fit$best <- 5}
	}

result$out$largestpatch = mean(lpatches)
result$out$largestpatch_sd = sd(lpatches)

result$out$best = c("DES","PL", "TPLup", "TPLdown", "EXP", "COV")[result$fit$best+1]
if(! result$out$best %in% c("DES","COV")) {
	result$out$p1 = result$fit$summary[[result$fit$best]][[2]]
	result$out$p2 = result$fit$summary[[result$fit$best]][[3]]
	result$out$p3 = result$fit$summary[[result$fit$best]][[4]]
	result$out$p1_sd = result$fit$summary[[result$fit$best]][[5]]
	result$out$p2_sd = result$fit$summary[[result$fit$best]][[6]]
	result$out$p3_sd = result$fit$summary[[result$fit$best]][[7]]
} else {
	result$out$p1 = NA
	result$out$p2 = NA
	result$out$p3 = NA
	result$out$p1_sd = NA
	result$out$p2_sd = NA
	result$out$p3_sd = NA
}

			# save result to file
	save(result, file = paste("results/", i, sep = ""))
	
#	Sys.sleep(.4)


}


filenames <- list.files("results/")
cuts <- round(seq(0,length(filenames), length = 12))

foreach(j = 1:11, .combine = "rbind") %dopar% {
files <- filenames[(cuts[j]+1):cuts[j+1]]

output <- data.frame(
			ID = NA,
			starting = NA, 
			globalgrazing = NA,
			stock = NA,
			g = NA, 
			b = NA, 
			m0 = NA,
			mortality = NA,
			mortality_border = NA,
			rho_plus = NA, 
			rho_plus_sd = NA,
			rho_plus_ini = NA,
			q_plus = NA,
			rho_zero = NA,
			rho_minus = NA,
			largestpatch = NA, 
			largestpatch_sd = NA,
			best =  NA,
			p1 = NA,
			p2 = NA,
			p3 = NA,
			p1_sd = NA,
			p2_sd = NA,
			p3_sd = NA,
			stable = NA
	)[-1,]


	
for(i in files) {

load(paste("results\\", i, sep = ""))

output <- rbind(output, result$out)
}

return(output)
} -> output
write.table(output, "output.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)



stopCluster(cl)
