#################################################
#
#		model: facilitation model plus grazing
#
#
#author: Flo
#date: 12.09.2013
#
#todo: funktionen auslagern
#  		gradienten definieren
#################################################

rm(list=ls())

########################################################################################
	

count  <- function(x, neighbor, state = NULL) {
			if(is.null(state)) state = neighbor
			
			neighbors <- numeric(length = prod(x$dim))
			x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
			for(k in interact) {
				neighbors <- neighbors + x_logical_with_border[x_to_evaluate+k]
			}
			return(neighbors)	
}

# function which counts the density of a state within neighboring cells of a particular state
# requires a previous global definition of the transformation vectors "x_to_evaluate" and "x_with_border" and the neighborhood vector "interact"
localdens <- function(x, neighbor, state) {
			neighbors <- numeric(length = length(which(x$cells == state))) # preallocate memory
			
			# create a logical vector x(which cells are in state "state"?) 
			# and transfer the vector x to a vector x_with_border: x[x_with_border]
			x_logical_with_border <- (x$cells == neighbor)[x_with_border]
			
			#of those values which are original values (x_to_evaluate) and which are "state": x_to_evaluate[initial$cells == "state"]
			x_state <- x_to_evaluate[x$cells == state]
			
			if(length(x_state) > 0) { 
			for(k in interact) {
				neighbors <- neighbors + x_logical_with_border[x_state+k]
			}
			# derive a vector which contains the value to count (also "state") in one of the neighbor cells j
			# cumulate the result for each of the four neighboring cells (for(k in interact))
			} else neighbors <- 0
			
			
			mean(neighbors/4, na.rm = TRUE)	
			# devide by 4 to get density
			# average over all cells
}


# get patch size and patchsize distribution
patches <- function(x, state, cumulative = TRUE) {
	pattern <- x$cells
	pattern <- pattern %in% state
	map <- rep(NA, times = prod(x$dim))
	old <- rep(99, times = prod(x$dim)) 
	
	while(!identical(old[pattern], map[pattern])) {
		old <- map
		count = as.integer(1)
		for(i in which(pattern)) {
			neighbors <- map[x_with_border][x_to_evaluate[i]+interact]
			if(all(is.na(neighbors)) ) { 
				map[i] <- count
			} else {
				map[i] <- min(neighbors, na.rm = TRUE)
			}
				count <- count +1
			}

		}
	
	map <- as.factor(map)
	patchvec <- as.vector(sapply(levels(map), function(i) length(which(map == i) ) )) 
	
	out <- vector()
	if(length(patchvec) > 0) out <- sort(patchvec) else out <- NA
	#out <- data.frame(size = unique(out), freq = sapply(unique(out), function(i) length(which(out >= i)) ))
	return(out)
	
	} 


################ parameter settings
first_ID = 27541  #
global <- c(rep(TRUE, (81+81+61+61+41)), rep(TRUE, (61+61+61+41)), rep(FALSE, (81+81+61+41)), rep(FALSE, (61+61+41)) )
stock <- c(rep(FALSE, (81+81+61+61+41)), rep(TRUE, (61+61+61+41)), rep(FALSE, (81+81+61+41)), rep(TRUE, (61+61+41)) )

graz <- c(
			rep(0, 81), rep(0.025, 81), rep(0.05, 61), rep(0.075, 61), rep(.1, 41),
			rep(0.025, 61), rep(0.05, 61), rep(0.075, 61), rep(.1, 41),
			rep(0.025, 81), rep(0.05, 81), rep(0.075, 61), rep(.1, 41),
			rep(0.025, 61), rep(0.05, 61), rep(0.075, 41)
		)
env <- 	c( seq(0.1,0.3,0.0025) , seq(0.2,0.4,0.0025),  seq(0.4,0.55,0.0025), seq(0.5,0.65,0.0025), seq(0.65,0.75,0.0025), 
			seq(0.3,0.45,0.0025) , seq(0.5,0.65,0.0025),  seq(0.65,0.8,0.0025), seq(0.8,0.9,0.0025), 
			seq(0.25,0.45,0.0025) , seq(0.45,0.65,0.0025),  seq(0.65,0.8,0.0025), seq(0.9,1,0.0025), 
			seq(0.4,0.55,0.0025) , seq(0.6,0.75,0.0025),  seq(0.85,0.95,0.0025)
		)
envdone <- 	seq(0.2,1,.01) #seq(0.20,0.98,.02)+.01

init <- c(0.250, 0.271 , 0.831, 0.900) #0.271 , 0.831,

lgraz <- length(graz)
lenv <- length(env)
linit <- length(init)
lglobal <- length(global)
lstock <- length(stock)
replication <- 1

parameters <- data.frame(
		global = rep(global, each = linit), 
		stock = rep(stock, each = linit), 
		m = 0.05, 
		g = rep(graz, each = linit), 
		b = rep(env, each = linit), 
		starting = rep(init, times = lgraz), 
		d = 0.1, c_ = 0.2, del = 0.1, r = 0.01, f = 0.9)

parameters <- parameters[!parameters$b %in% envdone,] 

parameters <- cbind(ID = 1:dim(parameters)[1]+27540, parameters )

parameters <- parameters[which(parameters$g == 0.1 & parameters$global == parameters$stock),]

# specify lattice
width = 100
height = 100

# initial cell states
states = c("+","0","-")
prob = c(9/10,.9/10,0.1/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
timesteps = 2000
delta = 1/5

t_eval <- ((timesteps/delta)-1000/delta):(timesteps/delta)+1
	
# derive helper vectors for counting: 
# transformation vector for evaluation at the border of the grid
	# set evaluation matrix 
	X <- matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
	# setting the border of the evaluation matrix X
	X <- cbind(X[,width], X, X[,1] )  
	X <- rbind(X[height,], X, X[1,] ) 
	# transformation vector which adds the border to the lattice:
	x_with_border <- as.integer(t(X))
	
	# from the matrix X (with copied border cells), which cells are the actual cells (reverse transformation vector of the previous lines) 
	x_to_evaluate <- sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]	)		

# defining the neighborhood which is to be evaluated	
	# set interaction matrix
	I <- matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)	
	# coordinates of neighbours in Interaction matrix I: 
	neighbours_in_I <- which(is.finite(abs(I)/abs(I)), arr.in = TRUE)
	# coordinates relative to the evaluated cell (=  which(is.na(I) ) 
	relrow <- neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
	relcol <- neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
	
	# relative position of the four direct neighbours of a cell
	interact <- (relrow * dim(X)[2] + relcol)

	
	
################ output saving

snapshots <-  (timesteps-20:0*50)/delta+1

#dir.create("results")

################ starting parallel backend
	
library(foreach)
library(doSNOW)

#workerlist <- rep(list(ubuWorker), times = 2)
workerlist <- rep("localhost", times = 15)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)


filenames <- list.files("results/")

################ starting foreach loop
foreach(iteration = parameters$ID) %dopar% {

#iteration = 30515 # 20063
set.seed(iteration)

# sampling the initial random grid into a list object
initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = prob ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)


# pre-run of 1000 timesteps without grazing (g = 0) and high or low environemntal parameter. 

 parms_temp <- list(
		m = 0.1, 		# intrinsic mortality
		b = parameters[parameters$ID == iteration,"starting"] , 		# beta*eps 
		d = 0.1,		# degradation
		c_ = 0.2, 		# beta*g  
		del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
		r = 0.01, 	# regeneration rate
		f = 0.9, 		# local fascilitation
		g = 0		#grazing
		)
	
for(i in 1:(500/delta+1)) {    #calculation loop
 if(i == 1) x_old <- initial
 parms_temp$rho <- sum(x_old$cells == "+")/(width*height)
 if(parms_temp$rho == 0) {flag <- FALSE} else {flag <- TRUE}
 
 x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

 if(flag) {	
		parms_temp$Q_plus <- count(x_old, "+")/4   	

		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		death <- with(parms_temp, (m*delta))
		death[death > 1] <- 1
		degradation <- with(parms_temp, (d *delta))
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in pre-run of", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in pre-run of", iteration, "in time step", i, "! balance parameters!!!")) 

		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"
	}  else {
		stop(paste("pre-run causes extinction in iteration", iteration, "! set higher environmental parameter b!"))
	}


		x_old <- x_new
		
} # end of pre-run simulation.

	x_0 <- x_new

# initialising result list 
result <- list()   # create an output list and
			result$time <- seq(0, timesteps, delta) # write simulated timesteps
			
			result$mortality <- numeric(length = timesteps/delta+1)
			result$mortality_border <- numeric(length = timesteps/delta+1)
			
			#save all level global densities: rho_*
			result$rho  <- list()			
			for(i in 1:length(states)) {
				result$rho[[i]] <- numeric(length = timesteps/delta+1)
				result$rho[[i]][1]  <- sum(x_0$cells == states[i])/(width*height) # write initial rho 
			}
			
			result$q_  <- list()			
			
			
			#for(j in 1:length(states)) {
			result$q_[[1]] <- numeric(length = timesteps/delta+1)
			result$q_[[1]][1]  <- localdens(x_0, "+", "+")# write initial local densities 
			
			#}

#  initialising first iteration 
	x_old <- x_0    #old landscape for first timestep   
	
		# copying parms for the simulation 
	#write.table(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, i$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
	result$timeseries <- list()
	result$timeseries[[1]] <- data.frame(
			ID = iteration, 
			global = TRUE, 
			stock = TRUE, 
			starting_b =  parms_temp$b, 
			timestep = 1, 
			g = parms_temp$g, 
			m0 = parms_temp$m , 
			b = parms_temp$b, 
			rho_plus = result$rho[[1]][1], 
			q_plus = result$q_[[1]][1])

	result$timeseries[[2]] <- x_0
	
	result$patches <- list()
	result$patches[[1]] <- patches(x_0,"+") 
	
	
	parms_temp <- as.list(parameters[parameters$ID == iteration,])

for(i in 2:(timesteps/delta+1)) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 
		
		if(parms_temp$rho == 0) {flag <- FALSE} else {flag <- TRUE}
		
		
	 # count local density of occupied fields for each cell: 
		parms_temp$Q_plus <- count(x_old, "+")/4   

# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
	if(flag) { # if density is unequal 0, then

		# calculate recolonisation rates of all cells
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		
		# calculate death rates
			# switch for livestock model. 
			rho_comp <- c(0.5, parms_temp$rho)[parms_temp$stock+1]
			
			# for differentiated grazing, i.e. associative protection, 
		if(parms_temp$global == FALSE) { 
				# determine vulnerability
				parms_temp$vul <- sum((1-parms_temp$Q_plus)[x_new$cells == "+" & (1-parms_temp$Q_plus) > 0]) /  sum(x_new$cells == "+")#/(width*height)

				death <- with(parms_temp, (m+g/rho_comp*(1-Q_plus)/vul)*delta)
			
		} else {
			# for undifferentiated grazing

				death <- with(parms_temp, (m+g/rho_comp)*delta)				
		}
		death[death > 1] <- 1
	
		}
		
		degradation <- with(parms_temp, (d *delta))
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		# check for sum of probabilities to be inferior 1 and superior 0
		if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		# apply rules 
	if(flag) {	
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		}
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		

		result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)/delta
		
		result$mortality_border[i] <- mean(death[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0], na.rm = TRUE)/delta
		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}


		result$q_[[1]][i]  <- localdens(x_new, states[j], "+")
		

		if(i %in% snapshots ) {
				result$timeseries[[1]] <- rbind(result$timeseries[[1]], data.frame(
										ID = iteration, 
										global = parms_temp$global, 
										stock = parms_temp$stock, 
										starting_b = parameters[parameters$ID == iteration,"starting"], 
										timestep = (i-1)*delta, 
										g = parms_temp$g, 
										m0 = parms_temp$m , 
										b = parms_temp$b, 
										rho_plus =  result$rho[[1]][i], 
										q_plus = result$q_[[1]][i])
								)

				result$timeseries[[length(result$timeseries)+1]] <- x_new 
				result$patches[[length(result$patches)+1]] <- patches(x_new,"+")
						
		}

		# activate to save timeseries at activation of herbivory
		#if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		#if(i %in% seq(1,timesteps, 50)/delta) write.table(c(iteration, i*delta, flag, x_new$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		x_old <- x_new
		
	} # end of simulation.
	
	
	
	##### log-binned patch-size distribution
	logbins <-  10^seq(log10(1),log10(10000), .1)
	result$logbin <- data.frame(size = logbins)
	
	for(j in 2:length(result$patches)) {
		result$logbin[paste(((c(NA,snapshots)-1)*delta)[j])] <- sapply(1:length(logbins), function(k) length(which(result$patches[[j]] >= logbins[k] & result$patches[[j]] < logbins[k+1])) )
	}
	
	result$logbin$mean <- apply(result$logbin[,2:length(result$patches)], 1, mean)
	result$logbin$sd <- apply(result$logbin[,2:length(result$patches)], 1, sd)
	
	
	#### cumulative patch size distribution
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


result$fit <- list()
result$fit$PL <- data.frame(snapshot = NA, a = NA, alpha = NA, third = NA, p_a = NA, p_alpha = NA, p_third = NA, AIC = NA)[-1,]
result$fit$TPLup <- data.frame(snapshot = NA, a = NA, alpha = NA, b = NA, p_a = NA, p_alpha = NA, p_b = NA, AIC = NA)[-1,]
result$fit$TPLdown <- data.frame(snapshot = NA, a = NA, alpha = NA, Sx = NA, p_a = NA, p_alpha = NA, p_Sx = NA, AIC = NA)[-1,]
result$fit$EXP <- data.frame(snapshot = NA, a = NA, eps = NA, third = NA, p_a = NA, p_alpha = NA, p_third = NA,  AIC = NA)[-1,]
result$fit$best <- vector()

	for(j in 1:length(snapshots)+1 ) {
	
	#j = 2
	
	if(is.na(result$cumpatch[[j]]) || dim(result$cumpatch[[j]])[1] < 3) {flag = FALSE} else { flag = TRUE }

if(flag) {
	if(!is.na(result$cumpatch[[j]][1,1])) dd3 <- data.frame(size = result$cumpatch[[j]]$size, n = result$cumpatch[[j]]$p)

	models  <- list()
	models$AIC <- vector()
	#############
PLlm <- lm(I(log(n)) ~  I(log(size)) , data = dd3) 
	
try({models$PL <- nls(I(log(n)) ~ log(a) - alpha * log(size), 
		data = dd3,
		start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]),
		trace = FALSE,
		nls.control(maxiter = 100)
		)}, silent = TRUE
	)
	
	
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
	

		
try( {models$TPLdown <- nls(I(log(n)) ~ I( log(a) - alpha * log(size) - (size/Sx) ), 
		data = dd3,
		start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2], Sx = 100),
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
	
	
	result$fit$best[j-1] <- which.min(models$AIC)

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
			rho_plus_ini = mean(x_0$cells == "+"),   #sum(x_0$cells == "+")/(width*height)		
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
	save(result, file = paste("results/result", iteration, sep = "_"))
	
#write.table(out, "output.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
#write.table(result$timeseries, "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)


gc() 
}


filenames <- list.files("results/")
cuts <- round(seq(0,length(filenames), length = 24))


foreach(j = 1:23, .combine = "rbind") %dopar% {
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

load(paste("results/", i, sep = ""))

output <- rbind(output, result$out)

}

return(output)
} -> output
write.table(output, "output.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)

foreach(j = 1:23, .combine = "rbind") %dopar% {
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
