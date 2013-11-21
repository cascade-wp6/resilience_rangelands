
#setwd("E:\\Eigene Dokumente\\Uni\\projects\\2013 CASCADE grazing\\sim9\\")
setwd("C:\\Users\\SCHNEIDER\\Documents\\projects\\CASCADE\\2013 CASCADE grazing\\sim9\\")



output <- read.csv("output.csv")
output <- rbind(output, read.csv("output2.csv"))



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
#i = 2360 
# output$ID[which(output$rho_plus > 0.02  & output$rho_plus < 0.1)] 
#i = 900# 328 # 2220 # 

load(paste("results\\result_", i, sep = ""))

	logbins <-  10^seq(log10(1),log10(10000), .1)
	result$logbin <- data.frame(size = logbins)
	
	for(j in 2:length(result$patches)) {
		result$logbin[paste(((c(NA,snapshots)-1)*delta)[j])] <- sapply(1:length(logbins), function(k) length(which(result$patches[[j]] >= logbins[k] & result$patches[[j]] < logbins[k+1])) )
	}
	
	result$logbin$mean <- apply(result$logbin[,2:length(result$patches)], 1, mean)
	result$logbin$sd <- apply(result$logbin[,2:length(result$patches)], 1, sd)
	
	result$cumpatch <- list()
	
for(j in 2:length(result$patches)) {

	if( !is.na(result$patch[j])) {
	cumbins <- sort(unique(unlist(result$patches[j]))) 
	#bins <- seq(1,10001, 10)

		result$cumpatch[[j]] <- data.frame(size = cumbins)
		result$cumpatch[[j]]$n <- sapply(cumbins, function(k) length(which(result$patches[[j]] >= k)) )
		result$cumpatch[[j]]$p <- sapply(cumbins, function(k) length(which(result$patches[[j]] >= k)) )/sum(result$cumpatch[[j]]$n)

	} else {
	result$cumpatch[[j]] <- NA
	}

}


#file.remove(paste("results\\", i, sep = ""))

result$fit <- list()
result$fit$PL <- data.frame(snapshot = NA, a = NA, alpha = NA, AIC = NA)[-1,]
result$fit$PLflat <- data.frame(snapshot = NA, a = NA, alpha = NA, b = NA, AIC = NA)[-1,]
result$fit$TPL <- data.frame(snapshot = NA, a = NA, alpha = NA, Sx = NA, AIC = NA)[-1,]
result$fit$EXP <- data.frame(snapshot = NA, a = NA, eps = NA,  AIC = NA)[-1,]
result$fit$best <- vector()

	for(j in 2:12 ) {
	
	#j = 2
	
	if(is.na(result$cumpatch[[j]]) || dim(result$cumpatch[[j]])[1] < 5) {flag = FALSE} else{ flag = TRUE }

if(flag) {
	
	dd3 <- data.frame(size = result$cumpatch[[j]]$size, n = result$cumpatch[[j]]$p)

	models  <- list()
	models$AIC <- vector()
	#############
models$PL <- lm(I(log(n)) ~  I(log(size)) , data = dd3) 
	
try({models$PLnls <- nls(I(log(n)) ~ log(a) - alpha * log(size), 
		data = dd3,
		start = list(a = exp(models$PL$coefficients[1]), alpha =  -models$PL$coefficients[2]),
		trace = FALSE,
		nls.control(maxiter = 100)
		)}, silent = TRUE
	)
	
	
	if(!is.null(models$PLnls)) {
		models$AIC[1] <- AIC(models$PLnls)
		result$fit$PL[j-1,] <- data.frame(snapshot = (snapshots[j-1]-1)*delta, a = coefficients(models$PLnls)[1], alpha = coefficients(models$PLnls)[2], AIC = AIC(models$PLnls))
	} else {
		models$PLnls  <- list(NA)
		models$AIC[1] <- NA
	}

###########
b=1/sum(result$cumpatch[[j]]$n)
try({models$PLflat <- nls(I(log(n)) ~ -I(log(a)) - alpha * log(size) + log(1+b/(a*size^(-alpha))), 
		data = dd3,
		start = list(a = 1, alpha =  4), #, b = 1/sum(result$cumpatch[[j]]$n)
        trace = FALSE,
		nls.control(maxiter = 50)
		)}, silent = TRUE
	)
	
	if(!is.null(models$PLflat)) {
		models$AIC[2] <- AIC(models$PLflat) 
		result$fit$PLflat[j-1,] <- data.frame(snapshot = (snapshots[j-1]-1)*delta, a = coefficients(models$PLflat)[1], alpha = coefficients(models$PLflat)[2], b = b, AIC = AIC(models$PLflat))

	} else { 
		models$PLflat  <- list(NA)
		models$AIC[2] <- NA
	}

###########
		
try( {models$TPLnls <- nls(I(log(n)) ~ -I(log(a)) - alpha * I(log(size)) - I(size/Sx), 
		data = dd3,
		start = list(a = exp(models$PL$coefficients[1]), alpha =  -models$PL$coefficients[2], Sx = 100),
        trace = FALSE
		)}, silent = TRUE
	)		

	if(!is.null(models$TPLnls)) {
		models$AIC[3] <- AIC(models$TPLnls) 
		result$fit$TPL[j-1,] <- data.frame(snapshot = (snapshots[j-1]-1)*delta, a = coefficients(models$TPLnls)[1], alpha = coefficients(models$TPLnls)[2], Sx = coefficients(models$TPLnls)[3], AIC = AIC(models$TPLnls))

	} else {
		models$TPLnls <- list(NA)
		models$AIC[3] <- NA
	}

###########
	
try( {models$EXP <- nls(I(log(n)) ~ -I(log(a)) -(eps*size) , 
		data = dd3,
		start = list(a = 1 ,eps = 0.8348),
        trace = FALSE
		)}, silent = TRUE
	)
	
		
	if(!is.null(models$EXP)) {
		models$AIC[4] <- AIC(models$EXP) 
		result$fit$EXP[j-1,] <- data.frame(snapshot = (snapshots[j-1]-1)*delta, a = coefficients(models$EXP)[1], eps = coefficients(models$EXP)[2], AIC = AIC(models$EXP))

	} else {
		models$EXP <- list(NA)
		models$AIC[4] <- NA
	}
	
	models$dAIC <- 	models$AIC -min(models$AIC, na.rm = TRUE)
	
	
	result$fit$best[j-1] <- which.min(models$AIC)

	} else {

	result$fit$best[j-1] <- NA

	}


}


best =  as.numeric(names(table(result$fit$best))[which.max(table(result$fit$best))])

if(length(best)) {
	result$fit$out <-  data.frame(
	ID = i,
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
	ID = i,
	best =  NA,
	p1 = NA,
	p2 = NA,
	p3 = NA,
	p1_sd = NA,
	p2_sd = NA,
	p3_sd = NA
	)
	
}
save(result, file = (paste("results\\result_", i, sep = "")) )

#	Sys.sleep(.4)


}

stopCluster(cl)


out_fit <- data.frame(
	ID = NA,
	best =  NA,
	p1 = NA,
	p2 = NA,
	p3 = NA,
	p1_sd = NA,
	p2_sd = NA,
	p3_sd = NA
	)[-1,]
	
for(i in output$ID) {

load(paste("results\\result_", i, sep = ""))

out_fit <- rbind(out_fit, result$fit$out)

}

write.table(out_fit, "output_fit.csv", row.names = FALSE, col.names = TRUE, sep = ",", append = TRUE)

output <- cbind(output[order(output$ID),] , out_fit[order(out_fit$ID),-1])

