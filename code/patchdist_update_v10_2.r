
#setwd("E:\\Eigene Dokumente\\Uni\\projects\\2013 CASCADE grazing\\sim9\\")
setwd("C:\\Users\\SCHNEIDER\\Documents\\projects\\CAS01_grazing\\data\\sim10\\")



# time and resolution of simulation
timesteps = 2000
delta = 1/5

snapshots <-  (timesteps-20:0*50)/delta+1

output <- read.csv("output.csv")



library(foreach)
library(doSNOW)

workerlist <- rep("localhost", times = 11)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)


foreach(i = output$ID) %dopar% {
#i =   2603  #which(! 1:2640  %in% read.csv("output_fit.csv")$ID)

# output$ID[which(output$rho_plus > 0.4  & output$rho_plus < 0.5)] 
#i = 1405
# output$ID[which(output$rho_plus > 0.68  & output$rho_plus < 0.7)] 
#i =  25508
# output$ID[which(output$rho_plus > 0.02  & output$rho_plus < 0.1)] 
#i = 900# 328 # 2220 # 
#i = 29880 #5053 #
#i = 4749
load(paste("results\\result_", i, sep = ""))

#dd4 <- do.call("rbind", result$cumpatch)
#file.remove(paste("results\\", i, sep = ""))



result$fit <- list()
result$fit$PL <- data.frame(snapshot = NA, a = NA, alpha = NA, third = NA, p_a = NA, p_alpha = NA, p_third = NA, AIC = NA)[-1,]
result$fit$TPLup <- data.frame(snapshot = NA, a = NA, alpha = NA, b = NA, p_a = NA, p_alpha = NA, p_b = NA, AIC = NA)[-1,]
result$fit$TPLdown <- data.frame(snapshot = NA, a = NA, alpha = NA, Sx = NA, p_a = NA, p_alpha = NA, p_Sx = NA, AIC = NA)[-1,]
result$fit$EXP <- data.frame(snapshot = NA, a = NA, eps = NA, third = NA, p_a = NA, p_alpha = NA, p_third = NA,  AIC = NA)[-1,]
result$fit$best <- vector()

	for(j in 2:22 ) {
	
	#j = 6
	
	if(is.na(result$cumpatch[[j]]) || dim(result$cumpatch[[j]])[1] < 3) {flag = FALSE} else { flag = TRUE }

if(flag) {
	if(!is.na(result$cumpatch[[j]][1,1])) dd3 <- data.frame(size = result$cumpatch[[j]]$size, n = result$cumpatch[[j]]$p)


	models  <- list()
	models$AIC <- vector()
	#############
PLlm <- lm(I(log(n)) ~  I(log(size)) , data = dd3) 

#if( max(dd3$size) <= 7000) {
try({models$PL <- nls(I(log(n)) ~ log(a) - alpha * log(size), 
		data = dd3,
		start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]),
		trace = FALSE,
		#algorithm = "port",
		nls.control(maxiter = 100)
		)}, silent = TRUE
	)
#	} 
	
	if(!is.null(models$PL)) {
		models$AIC[1] <- AIC(models$PL)
		result$fit$PL[j-1,] <- data.frame(snapshot = (snapshots[j-1]-1)*delta, a = coefficients(models$PL)[1], alpha = coefficients(models$PL)[2], third = NA, p_a = summary(models$PL)$coefficients[1,4], p_alpha = summary(models$PL)$coefficients[2,4], p_third = NA, AIC = AIC(models$PL))
	} else {
		models$PL  <- list(NA)
		models$AIC[1] <- NA
	}

###########

b=result$cumpatch[[j]]$p[dim(result$cumpatch[[j]])[1]] #1/sum(result$cumpatch[[j]]$n)  
try({models$TPLup <- nls(I(log(n)) ~ I( log(a) - alpha * log(size) + log(1+b/(a*size^(-alpha))) ), 
		data = dd3,
		start = list(a =  exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]), #, b = 1/sum(result$cumpatch[[j]]$n)
        trace = FALSE,
		#algorithm = "port",
		#lower = c(0, 0), upper = c(1, NA),
		nls.control(maxiter = 50)
		)}, silent = TRUE
	)
	
	
	if(!is.null(models$TPLup)) {
		models$AIC[2] <- AIC(models$TPLup) 
		result$fit$TPLup[j-1,] <- data.frame(snapshot = (snapshots[j-1]-1)*delta, a = coefficients(models$TPLup)[1], alpha = coefficients(models$TPLup)[2], b = b, p_a = summary(models$TPLup)$coefficients[1,4], p_alpha = summary(models$TPLup)$coefficients[2,4], p_b = NA, AIC = AIC(models$TPLup))

	} else { 
		models$TPLup  <- list(NA)
		models$AIC[2] <- NA
	}
	
	
try( {models$TPLdown <- nls(I(log(n)) ~ I( log(a) - alpha * log(size) - (size * Sx) ), 
		data = dd3,
		start = list(a = dd3$n[1], alpha =  exp(PLlm$coefficients[2]), Sx = 1/100),
        #algorithm = "port",
		trace = FALSE
		)}, silent = TRUE
	)		

	if(!is.null(models$TPLdown)) {
		models$AIC[3] <- AIC(models$TPLdown) 
		result$fit$TPLdown[j-1,] <- data.frame(snapshot = (snapshots[j-1]-1)*delta, a = coefficients(models$TPLdown)[1], alpha = coefficients(models$TPLdown)[2], Sx = coefficients(models$TPLdown)[3], p_a = summary(models$TPLdown)$coefficients[1,4], p_alpha = summary(models$TPLdown)$coefficients[2,4], p_Sx = summary(models$TPLdown)$coefficients[3,4], AIC = AIC(models$TPLdown))

	} else {
		models$TPLdown <- list(NA)
		models$AIC[3] <- NA
	}

###########
	
try( {models$EXP <- nls(I(log(n)) ~ I(log(a) -(eps*size)) , 
		data = dd3,
		start = list(a = exp(PLlm$coefficients[1]) ,eps = 1),
        #algorithm = "port",
		trace = FALSE
		)}, silent = TRUE
	)
	
		
	if(!is.null(models$EXP)) {
		models$AIC[4] <- AIC(models$EXP) 
		result$fit$EXP[j-1,] <- data.frame(snapshot = (snapshots[j-1]-1)*delta, a = coefficients(models$EXP)[1], eps = coefficients(models$EXP)[2], third = NA, p_a = summary(models$EXP)$coefficients[1,4], p_alpha = summary(models$EXP)$coefficients[2,4], p_Sx = NA, AIC = AIC(models$EXP))

	} else {
		models$EXP <- list(NA)
		models$AIC[4] <- NA
	}
	
		
	#min(models$AIC, na.rm = TRUE)
	models$dAIC <- 	models$AIC -min(models$AIC, na.rm = TRUE)
	
	
	result$fit$best[j-1] <- which.min(models$AIC[-4]+c(+0,0,0))

	} else {
	
		if(is.na(result$cumpatch[[j]])) { 
			result$fit$best[j-1] <- 0 
		} else { 
			result$fit$best[j-1] <- 5
		}

	}
}

result$fit$summary <- data.frame(
		ID = NA,
		snapshot = NA,	
		starting = NA, 
		globalgrazing = NA, 
		stock = NA,
		b = NA,
		g = NA,
		m0 = NA,
		largestpatch = NA,
		#rho_plus = ,
		bestmodel = NA,
		p1 = NA,
		p2 = NA,
		p3 = NA,		
		p1_p = NA,
		p2_p = NA,
		p3_p = NA
		)[-1,]

for(k in 1:length(result$fit$best) ) { 		# which(!result$fit$best %in% c(0,5))) {
result$fit$summary <- rbind(result$fit$summary, data.frame(
		ID = result$out$ID,
		snapshot = (snapshots[k]-1)*delta,	
		starting = result$out$starting, 
		globalgrazing = result$out$globalgrazing, 
		stock = result$out$stock,
		b = result$out$b,
		g = result$out$g,
		m0 = result$out$m0,
		largestpatch = max(result$patches[[k+1]]),
		#rho_plus = ,
		bestmodel = c("DES","PL", "TPLup", "TPLdown", "EXP", "COV")[result$fit$best[k]+1],
		p1 = NA,
		p2 = NA,
		p3 = NA,		
		p1_p = NA,
		p2_p = NA,
		p3_p = NA
		)
		)
		
	}

for(k in which(!result$fit$best %in% c(0,5))) {
	result$fit$summary[k,11:16] <- result$fit[[result$fit$best[k]]][k,2:7]
	}

best =  c("DES","PL", "TPLup", "TPLdown", "EXP", "COV")[as.numeric(names(table(result$fit$best))[which.max(table(result$fit$best))])+1]

	
if(! best %in% c("DES","COV")) {
	result$fit$out <-  data.frame(
	ID = result$out$ID,
	best =  best,
	p1 = mean(result$fit[[best]][,2], na.rm = TRUE),
	p2 = mean(result$fit[[best]][,3], na.rm = TRUE),
	p3 = mean(result$fit[[best]][,4], na.rm = TRUE),
	p1_sd = sd(result$fit[[best]][,2], na.rm = TRUE),
	p2_sd = sd(result$fit[[best]][,3], na.rm = TRUE),
	p3_sd = sd(result$fit[[best]][,4], na.rm = TRUE)
	)
} else {
	result$fit$out <-  data.frame(
	ID = result$out$ID,
	best =  best,
	p1 = NA,
	p2 = NA,
	p3 = NA,
	p1_sd = NA,
	p2_sd = NA,
	p3_sd = NA
	)
	
}


result$out$largestpatch = res
result$out$largestpatch_sd = 
result$out$best = result$fit$out$best
result$out$p1 = result$fit$out$p1
result$out$p2 = result$fit$out$p2
result$out$p3 = result$fit$out$p3
result$out$p1_sd = result$fit$out$p1_sd
result$out$p2_sd = result$fit$out$p2_sd
result$out$p3_sd = result$fit$out$p3_sd

			# save result to file
	save(result, file = paste("results/result", result$out$ID, sep = "_"))
	
#	Sys.sleep(.4)


}


filenames <- list.files("results/")
cuts <- round(seq(0,length(filenames), length = 12))

foreach(j = 1:11, .combine = "rbind") %dopar% {
files <- filenames[(cuts[j]+1):cuts[j+1]]

out_fits <- data.frame(
		ID = NA,
		snapshot = NA,	
		starting = NA, 
		globalgrazing = NA, 
		stock = NA,
		b = NA,
		g = NA,
		m0 = NA,
		largestpatch = NA,
		bestmodel = NA,
		p1 = NA,
		p2 = NA,
		p3 = NA,		
		p1_p = NA,
		p2_p = NA,
		p3_p = NA
	)[-1,]

	
for(i in files) {

load(paste("results/", i, sep = ""))

out_fits <- rbind(out_fits, result$fit$summary)

}

return(out_fits)
} -> out_fits
write.table(out_fits, "output_fits.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)

stopCluster(cl)
