#################################################
#
#		model: facilitation model (Kefi et al 2007)
#
#
#author: Flo
#date: 02.04.2013
#
#
#	Tasks for Flo: 
#		identify bottlenecks for speedy and parallel computing
#		implement test for stable equilibrium
#		
#   Tasks for Marina:
# 		adapt for two-species model (implement rules, find parameters)
# 		adapt output
#
#################################################

rm(list=ls())

setwd("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\thread herbivory\\")

env <- 	seq(0.2,1,.05)
graz <- seq(0,.9, 0.3)
mort <- seq(0.1,.3, 0.1)
init <- c("low", "high")

lgraz <- length(graz)
lenv <- length(env)
lmort <- length(mort)

env <-	rep(rep(env, times = lgraz), times = lmort)
graz <-	 rep(rep(graz, each = lenv), times = lmort)
mort <-  rep(mort, each = lgraz*lenv)

replicates = 2
init <- rep(init, each = length(mort) )
env <- rep(env, times = replicates)
graz <- rep(graz, times = replicates)
mort <- rep(mort, times = replicates)

# defining parameter set
parameters = data.frame(
	ID = 1:length(graz), 
	starting = init,
	g = graz,		#grazing
	m = mort, 		# intrinsic mortality
	b = env, 		# beta*eps 
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.01, 	# regeneration rate
	f = 0.9		# local fascilitation
)


# defining parameter set
parameters_ini = list(
	low = list(
		m = 0.1, 		# intrinsic mortality
		b = 0.25, 		# beta*eps 
		d = 0.1,		# degradation
		c_ = 0.2, 		# beta*g  
		del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
		r = 0.01, 	# regeneration rate
		f = 0.9, 		# local fascilitation
		g = 0		#grazing
		), 
	high = list(
		m = 0.1, 		# intrinsic mortality
		b = 0.9, 		# beta*eps 
		d = 0.1,		# degradation
		c_ = 0.2, 		# beta*g  
		del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
		r = 0.01, 	# regeneration rate
		f = 0.9, 		# local fascilitation
		g = 0		#grazing
	)
	)


# specify lattice
width = 10
height = 10

# initial cell states
states = c("+","0","-")
prob = c(9/10,.9/10,0.1/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
timesteps = 200
addgrazing = 50
delta = 1/10



out_header <- data.frame(
			ID = NA,
			starting = NA, 
			grazing = NA, 
			environmt = NA, 
			mort_int = NA, 			
			mortality = NA,
			mortality_border = NA,
			rho_plus = NA, 
			rho_plus_ungrazed = NA, 
			q_plus = NA,
			q_plus_ungrazed = NA,
			rho_zero = NA,
			rho_minus = NA
			)

write.table(out_header[-1,], "output.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)



header_grids <- c("ID", "grazing", "starting", "timestep", "g", "m" , "b", paste("X", 1:(width*height), sep =""))
snapshots <- c(1, c(addgrazing, timesteps-4:1*50, timesteps)/delta+1)



count  <- function(x, neighbor, state = NULL) {
			if(is.null(state)) state = neighbor
			
			neighbors <- numeric(length = prod(x$dim))
			x_logical_with_border <- (x$cells == neighbor)[x_with_border]
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

########################################################################################
########################################################################################
	
	
library(foreach)
library(doSNOW)

#	, rep(list(ubuWorker), times = 23)
workerlist <- rep("localhost", times = 11)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)


foreach(iteration = 1:length(env)) %dopar% {


#### first run with grazing
# sampling the initial random grid into a list object
initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = prob ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)


# initialising result list 
result <- list()   # create an output list and
			result$time <- seq(0, timesteps, delta) # write simulated timesteps
			
			result$mortality <- numeric(length = timesteps/delta+1)
			#save all level global densities: rho_*
			result$rho  <- list()			
			for(i in 1:length(states)) {
				result$rho[[i]] <- numeric(length = timesteps/delta+1)
				result$rho[[i]][1]  <- sum(initial$cells == states[i])/(width*height) # write initial rho 
			}
			
			result$q_  <- list()			
			for(j in 1:length(states)) {
			result$q_[[j]] <- numeric(length = timesteps/delta)
			result$q_[[j]][1]  <- mean(rowSums( sapply(interact, function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == "+"]+k]) )) /4# write initial local densities 
			}
			


#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	
	parms_temp <- as.list(parameters[iteration,])
		# copying parms for the simulation 
	#write.table(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, i$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
	if(1 %in% snapshots) write.table(as.data.frame(matrix(c( parms_temp$ID, TRUE, parms_temp$starting, 1, parms_temp$g, parms_temp$m, parms_temp$b, initial$cells), nrow = 1, dimnames = list("1", header_grids))), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

	parms_temp <- parameters_ini[[parameters[iteration,2]]]

	
for(i in 2:(timesteps/delta+1)) {    #calculation loop

		if(i == addgrazing/delta ) parms_temp <- as.list(parameters[iteration,])

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 

	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		# This is a vectorised evaluation: first transform the original grid x_old, into a grid with border. Here check for all cells of the original grid (x_to_evaluate) the content of all four neighbouring cells, whose relative position is defined in interact. The number of cells with content "+" are count and divided by 4.
	# count local density of degraded fields for each cell
		parms_temp$Q_unveg <- ( count(x_old, c("0", "-")))/4 
		#parms_temp$vul <- sum(parms_temp$Q_unveg[x_old$cells == "+"]*4)/(width*height)
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		if(i < addgrazing/delta) {death <- with(parms_temp, m*delta)} else {death <- with(parms_temp, (m+g*Q_unveg)*delta)}
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > .95 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		

		result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)/delta
		
		result$mortality_border[i] <- mean(death[x_old$cells == "+" & parms_temp$Q_unveg > 0], na.rm = TRUE)/delta
		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}

		
		for(j in 1:length(states)) {
		result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		}
		
		if(i %in% snapshots) write.table(as.data.frame(matrix(c( parms_temp$ID, TRUE, parms_temp$starting, i, parms_temp$g, parms_temp$m, parms_temp$b, x_new$cells), nrow = 1, dimnames = list("1", header_grids))), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		# activate to save timeseries at activation of herbivory
		#if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		#if(i %in% seq(1,timesteps, 50)/delta) write.table(c(iteration, i*delta, flag, x_new$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		x_old <- x_new
		
		gc()  #garbage collection
	} # end of simulation.

	
	t_eval <- ((timesteps/delta)-50/delta):(timesteps/delta)+1
	
	out <- data.frame(
			ID = iteration,
			starting = parms_temp$starting, 
			grazing = parms_temp$g, 
			environmt = parms_temp$b, 
			mort_int = parms_temp$m, 			
			mortality = mean(result$mortality[t_eval], na.rm = TRUE),
			mortality_border = mean(result$mortality_border[t_eval], na.rm = TRUE),
			rho_plus = mean(result$rho[[1]][t_eval], na.rm = TRUE), 
			rho_plus_ungrazed = mean(result$rho[[1]][(t_eval-1)-(timesteps-addgrazing)/delta], na.rm = TRUE), 
			q_plus = mean(result$q_[[1]][t_eval], na.rm = TRUE),
			q_plus_ungrazed = mean(result$q_[[1]][(t_eval-1)-(timesteps-addgrazing)/delta], na.rm = TRUE),
			rho_zero = mean(result$rho[[2]][t_eval], na.rm = TRUE),
			rho_minus = mean(result$rho[[3]][t_eval], na.rm = TRUE)
			)

write.table(out, "output.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)






################################################################################ second run without grazing
# sampling the initial random grid into a list object
initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = prob ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)


# initialising result list 
result <- list()   # create an output list and
			result$time <- seq(0, timesteps, delta) # write simulated timesteps
			
			result$mortality <- numeric(length = timesteps/delta+1)
			#save all level global densities: rho_*
			result$rho  <- list()			
			for(i in 1:length(states)) {
				result$rho[[i]] <- numeric(length = timesteps/delta+1)
				result$rho[[i]][1]  <- sum(initial$cells == states[i])/(width*height) # write initial rho 
			}
			
			result$q_  <- list()			
			for(j in 1:length(states)) {
			result$q_[[j]] <- numeric(length = timesteps/delta)
			result$q_[[j]][1]  <- mean(rowSums( sapply(interact, function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == "+"]+k]) )) /4# write initial local densities 
			}
			


#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	
	parms_temp <- as.list(parameters[iteration,])
		# copying parms for the simulation 
	#write.table(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, i$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
	if(1 %in% snapshots) write.table(as.data.frame(matrix(c( parms_temp$ID, FALSE, parms_temp$starting, 1, 0, out$grazing, parms_temp$b, initial$cells), nrow = 1, dimnames = list("1", header_grids))), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

	parms_temp <- parameters_ini[[parameters[iteration,2]]]

	
for(i in 2:(timesteps/delta+1)) {    #calculation loop

		if(i == addgrazing/delta  ) {
			parms_temp <- as.list(parameters[iteration,])
			parms_temp$m <- out$grazing
			parms_temp$g <- 0
		}

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 

	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		# This is a vectorised evaluation: first transform the original grid x_old, into a grid with border. Here check for all cells of the original grid (x_to_evaluate) the content of all four neighbouring cells, whose relative position is defined in interact. The number of cells with content "+" are count and divided by 4.
	# count local density of degraded fields for each cell
		parms_temp$Q_unveg <- ( count(x_old, c("0", "-")))/4 
		#parms_temp$vul <- sum(parms_temp$Q_unveg[x_old$cells == "+"]*4)/(width*height)
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		if(i < addgrazing/delta) {death <- with(parms_temp, m*delta)} else {death <- with(parms_temp, (m+g*Q_unveg)*delta)}
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > .95 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		

		result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)/delta
		
		result$mortality_border[i] <- mean(death[x_old$cells == "+" & parms_temp$Q_unveg > 0], na.rm = TRUE)/delta
		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}

		
		for(j in 1:length(states)) {
		result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		}
		
		if(i %in% snapshots) write.table(as.data.frame(matrix(c( parms_temp$ID, TRUE, parms_temp$starting, i, parms_temp$g, parms_temp$m, parms_temp$b, x_new$cells), nrow = 1, dimnames = list("1", header_grids))), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		# activate to save timeseries at activation of herbivory
		#if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		#if(i %in% seq(1,timesteps, 50)/delta) write.table(c(iteration, i*delta, flag, x_new$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		x_old <- x_new
		
		gc()  #garbage collection
	} # end of simulation.

	
	t_eval <- ((timesteps/delta)-50/delta):(timesteps/delta)+1
	
	out <- data.frame(
			ID = iteration,
			starting = parms_temp$starting, 
			grazing = parms_temp$g, 
			environmt = parms_temp$b, 
			mort_int = parms_temp$m, 			
			mortality = mean(result$mortality[t_eval], na.rm = TRUE),
			mortality_border = mean(result$mortality_border[t_eval], na.rm = TRUE),
			rho_plus = mean(result$rho[[1]][t_eval], na.rm = TRUE), 
			rho_plus_ungrazed = mean(result$rho[[1]][(t_eval-1)-(timesteps-addgrazing)/delta], na.rm = TRUE), 
			q_plus = mean(result$q_[[1]][t_eval], na.rm = TRUE),
			q_plus_ungrazed = mean(result$q_[[1]][(t_eval-1)-(timesteps-addgrazing)/delta], na.rm = TRUE),
			rho_zero = mean(result$rho[[2]][t_eval], na.rm = TRUE),
			rho_minus = mean(result$rho[[3]][t_eval], na.rm = TRUE)
			)

write.table(out, "output.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)


gc() 
}


stopCluster(cl)
