# repair run 10


rm(list=ls())

##########


################ starting parallel backend
	
library(foreach)
library(doSNOW)

#workerlist <- rep(list(ubuWorker), times = 2)
workerlist <- rep("localhost", times = 23)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)

########################


################ parameter settings
first_ID = 1  #
env <- 	seq(0.2,1,.01) #seq(0.20,0.98,.02)+.01
graz <- seq(.0, .1, 0.025) #.025 #
mort <- 0.05 #c(0.05, 0.1) #seq(0.05,.15, 0.05)
init <- round(exp(seq(log(0.25), log(0.9), length = 33))[seq(1,33,2)], digits = 3)
global <- c(FALSE, TRUE)
stock <- c(FALSE, TRUE)


lgraz <- length(graz)
lenv <- length(env)
lmort <- length(mort)
linit <- length(init)
lglobal <- length(global)
lstock <- length(stock)
replication <- 1


# defining parameter set
parameters = data.frame(
	ID = first_ID:(lgraz*lmort*lglobal*linit*lenv*lstock*replication+first_ID-1), 
	global = rep(global, each = lgraz*lenv*lmort*linit*lstock*replication),
	stock =  rep(rep(stock, each = lgraz*lenv*lmort*linit*replication), times = lstock),
	starting = rep(rep(init, each = lgraz*lenv*lmort*replication), times = lglobal*lstock),
	m = rep(rep(mort, each = lgraz*lenv*replication), times = lglobal*linit*lstock), 		# intrinsic mortality
	g = rep(rep(graz, each = lenv*replication), times = lmort*lglobal*linit*lstock),		#grazing
	b = rep(rep(env, each = replication) , times = lgraz*lmort*lglobal*linit*lstock), 		# beta*eps 
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.01, 	# regeneration rate
	f = 0.9		# local fascilitation
)

# time and resolution of simulation
timesteps = 2000
delta = 1/5

t_eval <- ((timesteps/delta)-1000/delta):(timesteps/delta)+1
	

snapshots <-  (timesteps-20:0*50)/delta+1


filenames <- list.files("results/")

################ starting foreach loop
foreach(i = filenames) %dopar% {
#i = 1

load(paste("results/", i, sep = ""))
parms_temp <- as.list(parameters[parameters$ID == result$out$ID,])
iteration = result$out$ID

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
		ID = iteration,
		snapshot = (snapshots[k]-1)*delta,	
		starting = parms_temp$starting, 
		globalgrazing = parms_temp$global, 
		stock = parms_temp$stock,
		b = parms_temp$b,
		g = parms_temp$g,
		m0 = parms_temp$m,
		largestpatch = max(result$patches[[k]]),
		bestmodel = c("0","PL", "TPLup", "TPL", "EXP", "1")[result$fit$best[k]+1],
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

best =  as.numeric(names(table(result$fit$best))[which.max(table(result$fit$best))])

	
if(! best %in% c(0,5)) {
	result$fit$out <-  data.frame(
	ID = iteration,
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
	ID = iteration,
	best =  best,
	p1 = NA,
	p2 = NA,
	p3 = NA,
	p1_sd = NA,
	p2_sd = NA,
	p3_sd = NA
	)
	
}


	result$out <- data.frame(
			ID = iteration,
			starting = parms_temp$starting, 
			globalgrazing = parms_temp$global, 
			stock = parms_temp$stock,
			g = parms_temp$g,
			b = parms_temp$b, 
			m0 = parms_temp$m, 			
			mortality = mean(result$mortality[t_eval], na.rm = TRUE),
			mortality_border = mean(result$mortality_border[t_eval], na.rm = TRUE),
			rho_plus = mean(result$rho[[1]][t_eval], na.rm = TRUE), 
			rho_plus_sd = sd(result$rho[[1]][t_eval], na.rm = TRUE), 
			rho_plus_ini = mean(result$timeseries[[2]]$cells == "+"),   #sum(x_0$cells == "+")/(width*height)		
			q_plus = mean(result$q_[[1]][t_eval], na.rm = TRUE),
			rho_zero = mean(result$rho[[2]][t_eval], na.rm = TRUE),
			rho_minus = mean(result$rho[[3]][t_eval], na.rm = TRUE),
			largestpatch = mean(result$fit$summary$largestpatch, na.rm = TRUE),
			largestpatch_sd = sd(result$fit$summary$largestpatch, na.rm = TRUE),
			best = result$fit$out$best,
			p1 = result$fit$out$p1,
			p2 = result$fit$out$p2,
			p3 = result$fit$out$p3,
			p1_sd = result$fit$out$p1_sd,
			p2_sd = result$fit$out$p2_sd,
			p3_sd = result$fit$out$p3_sd,
			#largestpatch = , 
			stable = mean(result$rho[[1]][t_eval[1:((length(t_eval)-1)/2)] ]  ) - mean(result$rho[[1]][t_eval[((length(t_eval)+1)/2):length(t_eval)] ] ) < 0.01
	)
	# save result to file

save(result, file = paste("results/", i, sep = ""))


}



for(i in filenames) {

load(paste("results/", i, sep = ""))

output <- rbind(output, result$out)

out_fits <- rbind(out_fits, result$fit$summary)

}

write.table(output, "output.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)
write.table(out_fits, "output_fits.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)


stopCluster(cl)
