
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

# plotting function for objects of class "landscape".
plot.landscape <- function(x, grid = FALSE, axis = FALSE, cols = "auto", add = FALSE, ani = FALSE, ...) {
  lvls <- levels(x$cells) 
  nlev <- length(lvls)
  if(cols[1] == "auto") cols = c("black","grey80", "white")  # default color value
  
  if(ani & Sys.info()[['sysname']] == "Windows") adj = -0.5 else adj = 0 #this adjustment constant is added when producing a pixel accurate png or gif file, requires TRUE when the function is used to plot animated figures. 
  
  if(!add) plot(NA,NA, xlim = c(0.5+adj, x$dim[1]+0.5+adj), ylim = c( x$dim[2]+0.5+adj, 0+0.5+adj), bty = c("n", "o")[grid+1], xaxs = "i", yaxs = "i",xlab = "", ylab = "", xaxt = "n", yaxt = "n", ... ) 
  
  if(axis && !add) axis(3) 
  if(axis && !add) axis(2)
  
  if(grid) border = "grey80" else border = cols[as.numeric(x$cells)]
  
  rect(rep(1:x$dim[1], times = x$dim[2])-.5, rep(1:x$dim[2], each = x$dim[1])-.5, rep(1:x$dim[1], times = x$dim[2])+.5, rep(1:x$dim[2], each = x$dim[1])+.5, col = cols[as.numeric(x$cells)], border = border)
  
  if(grid) box()
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
first_ID = 1  #

global <- list(	
	m0 = 0.05, #intrinsic mortality
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.01, 	# regeneration rate
	f = 0.9		# local fascilitation
	)
	
parameters <- list(
	run = 1,
	b = seq(0.2,1,.001), #seq(0.20,0.98,.02)+.01
	g = 0, #.025 #
	assoc = c(FALSE),
	stock = c(FALSE),
	starting = seq(0.1,0.9, 0.4),
	iter = c(1)
)

iterations <- expand.grid(parameters)

parameters2 <- list(
	run = 2,
	b = seq(0.2,1,.002), #seq(0.20,0.98,.02)+.01
	g = c(0.025, 0.05, 0.075, 0.1), #.025 #
	assoc = c(FALSE, TRUE),
	stock = c(FALSE, TRUE),
	starting = seq(0.1,0.9, 0.4),
	iter = c(1)
)

iterations <- rbind(iterations, expand.grid(parameters2))

parameters3 <- list(
	run = 3,
	b = seq(0.2,1,.002), #seq(0.20,0.98,.02)+.01
	g = c(0.025, 0.05, 0.075, 0.1), #.025 #
	assoc = c(FALSE, TRUE),
	stock = c(FALSE, TRUE),
	starting = seq(0.1,0.9, 0.4),
	iter = c(1)
)

iterations <- rbind(iterations, expand.grid(parameters3))

iterations <- cbind(ID = 1:dim(iterations)[1],iterations, global)
str(iterations)



# specify lattice
width = 100
height = 100

# initial cell states
states = c("+","0","-")
prob = c(9/10,.9/10,0.1/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
timesteps = 3000
delta = 1/2

t_eval <- ((timesteps/delta)-1000/delta):(timesteps/delta)+1
	
snapshots <-  (timesteps-20:0*50)/delta+1

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

dir.create("results")

################ starting parallel backend
	
library(foreach)
library(doSNOW)

#workerlist <- rep(list(ubuWorker), times = 2)
workerlist <- rep("localhost", times = 23)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)


## first run sim for g = 0 in high res. with initial prob  drawn from uniform distribution. 
# output rho_plus_final
# collect database for each rho_plus levels with res 0.01
# define b gradient unambiguously with steplength 0.005 ? higher res if rho_plus is "spacy"?
# initialise all other sims with webs drawn randomly from the database. 
# run simulations over 3000 timesteps
# collect output
# check different levels of stability <0.01, <0.001 , <0.0001
# fit models



################ starting foreach loop
#foreach(iteration = iterations$ID[iterations$run == 1]) %dopar% {

foreach(iteration = c(1791, 943), .combine = "c") %do% {

#iteration =  1414
set.seed(iteration)

	parms_temp <- as.list(iterations[iterations$ID == iteration,])

# sampling the initial random grid into a list object
parms_temp$rho_plus <- abs(rnorm(1, mean = parms_temp$starting, sd = 0.2))
if(parms_temp$rho_plus > 1) parms_temp$rho_plus <- 1 - (parms_temp$rho_plus-1) 

initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = c(parms_temp$rho_plus, (1-parms_temp$rho_plus)/2, (1-parms_temp$rho_plus)/2 ) ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)
	
parms_temp$q_plus_plus <-  mean(count(initial, "+", "+")/4)

	result <- list()
	
	result$time <- seq(0, timesteps, delta) # write simulated timesteps
			
		#result$mortality <- numeric(length = timesteps/delta+1)
		#result$mortality_border <- numeric(length = timesteps/delta+1)
			
		#save all level global densities: rho_*
	result$rho_plus  <- vector("numeric", length = length(result$time))
			#for(i in 1:length(states)) {
			#	result$rho[[i]] <- numeric(length = timesteps/delta+1)
			#	result$rho[[i]][1]  <- sum(x_0$cells == states[i])/(width*height) # write initial rho 
			#}
			
		result$rho_plus[1] <- parms_temp$rho_plus

	result$q_plus_plus  <- vector("numeric", length = length(result$time))		
		result$q_plus_plus[1] <- parms_temp$q_plus_plus
						
			#for(j in 1:length(states)) {
			#result$q_[[1]] <- numeric(length = timesteps/delta+1)
			#result$q_[[1]][1]  <- localdens(x_0, "+", "+")# write initial local densities 
			#}

	result$snapshots <- data.frame(
			ID = parms_temp$ID, 
			assoc = parms_temp$assoc, 
			stock = parms_temp$stock, 
			timestep = 1, 
			g = parms_temp$g, 
			m0 = parms_temp$m0 , 
			b = parms_temp$b, 
			rho_plus = parms_temp$rho_plus, 
			q_plus_plus = parms_temp$q_plus_plus
			)

	result$grids <- list()
	result$grids[[1]] <- initial
	
	result$patches <- list()
	result$patches[[1]] <- patches(initial,"+") 
	
	x_old <- initial
	
for(i in 2:(timesteps/delta+1)) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		if(parms_temp$rho_plus == 0) {flag <- FALSE} else {flag <- TRUE}
		
	 # count local density of occupied fields for each cell: 
		parms_temp$Q_plus <- count(x_old, "+", "+")/4

# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
	if(flag) { # if density is unequal 0, then

		# calculate recolonisation rates of all cells
		recolonisation <- with(parms_temp, (del*rho_plus+(1-del)*Q_plus)*(b-c_*rho_plus)*delta)
		
		# calculate death rates
			# switch for livestock model. 
			rho_comp <- c(0.5, parms_temp$rho_plus)[parms_temp$stock+1]
			
			# for differentiated grazing, i.e. associative protection, 
		if(parms_temp$assoc == TRUE) { 
				# determine vulnerability
				parms_temp$vul <- sum((1-parms_temp$Q_plus)[x_new$cells == "+" & (1-parms_temp$Q_plus) > 0]) /  sum(x_new$cells == "+")#/(width*height)

				death <- with(parms_temp, (m0+g/rho_comp*(1-Q_plus)/vul)*delta)
			
		} else {
			# for undifferentiated grazing

				death <- with(parms_temp, (m0+g/rho_comp)*delta)				
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

		parms_temp$rho_plus <- sum(x_new$cells == "+")/(width*height) 
		result$rho_plus[i] <- parms_temp$rho_plus
		#result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)/delta
		#result$mortality_border[i] <- mean(death[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0], na.rm = TRUE)/delta
		
		parms_temp$rho_plus_plus <- mean(parms_temp$Q_plus[x_new$cells == "+"]/4) 
		result$rho_plus_plus[i] <- parms_temp$rho_plus_plus
		#for(j in 1:length(states)) {
		#	result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
		#	}


		#result$q_[[1]][i]  <- localdens(x_new, states[j], "+")
		

		if(i %in% snapshots ) {
				result$snapshots <- rbind(result$snapshots, data.frame(
										ID = iteration, 
										assoc = parms_temp$assoc, 
										stock = parms_temp$stock, 
										timestep = (i-1)*delta, 
										g = parms_temp$g, 
										m0 = parms_temp$m , 
										b = parms_temp$b, 
										rho_plus =  parms_temp$rho_plus, 
										q_plus_plus = parms_temp$rho_plus_plus
										)
								)

				result$grids[[length(result$grids)+1]] <- x_new 
				result$patches[[length(result$patches)+1]] <- patches(x_new,"+")
						
		}

		# activate to save timeseries at activation of herbivory
		#if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		#if(i %in% seq(1,timesteps, 50)/delta) write.table(c(iteration, i*delta, flag, x_new$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		x_old <- x_new
		
	} # end of simulation.

	result$out <- data.frame(
			ID = iteration,
			starting = parms_temp$starting, 
			assoc = parms_temp$assoc, 
			stock = parms_temp$stock,
			g = parms_temp$g,
			b = parms_temp$b, 
			m0 = parms_temp$m, 			
			#mortality = mean(result$mortality[t_eval], na.rm = TRUE),
			#mortality_border = mean(result$mortality_border[t_eval], na.rm = TRUE),
			rho_plus = mean(result$rho_plus[t_eval], na.rm = TRUE), 
			rho_plus_sd = sd(result$rho_plus[t_eval], na.rm = TRUE), 
			rho_plus_ini = mean(initial$cells == "+"),   #sum(x_0$cells == "+")/(width*height)		
			rho_plus_fin = parms_temp$rho_plus,   #sum(x_0$cells == "+")/(width*height)		
			stable1 = mean(result$rho_plus[t_eval[1:((length(t_eval)-1)/2)] ]  ) - mean(result$rho_plus[t_eval[((length(t_eval)+1)/2):length(t_eval)] ] ) < 0.001,
			stable2 = mean(result$rho_plus[t_eval[1:((length(t_eval)-1)/2)] ]  ) - mean(result$rho_plus[t_eval[((length(t_eval)+1)/2):length(t_eval)] ] ) < 0.01
	)
	
	save(result, file = paste("results/result", iteration, sep = "_"))

	return(result$grids[15:22])
} -> gridlist

# this is storing all the snaphot grids generated from the first run of the simulation
save(gridlist, file = "gridlist")

rho_plus_list <- sapply(gridlist, function(x) sum(x$cells == "+")/(width*height) )

initialgrids <- lapply(1:100, function(x) try(sample(which(rho_plus_list > seq(0,0.99,0.01)[x] & rho_plus_list < seq(0.01,1,0.01)[x]), 25, replace = TRUE) , silent = TRUE))






################ starting foreach loop
foreach(iteration = iterations$ID[iterations$run != 1]) %dopar% {

#iteration =  1414
set.seed(iteration)

	parms_temp <- as.list(iterations[iterations$ID == iteration,])

# sampling the initial random grid into a list object
parms_temp$rho_plus <- abs(rnorm(1, mean = parms_temp$starting, sd = 0.2))
if(parms_temp$rho_plus > 1) parms_temp$rho_plus <- 1 - (parms_temp$rho_plus-1) 

bin <- match(parms_temp$rho, seq(0,1,0.01))
initial <- gridlist[sample(initialgrids[[bin]], 1))]

parms_temp$q_plus_plus <-  mean(count(initial, "+", "+")/4)

	result <- list()
	
	result$time <- seq(0, timesteps, delta) # write simulated timesteps
			
		#result$mortality <- numeric(length = timesteps/delta+1)
		#result$mortality_border <- numeric(length = timesteps/delta+1)
			
		#save all level global densities: rho_*
	result$rho_plus  <- vector("numeric", length = length(result$time))
			#for(i in 1:length(states)) {
			#	result$rho[[i]] <- numeric(length = timesteps/delta+1)
			#	result$rho[[i]][1]  <- sum(x_0$cells == states[i])/(width*height) # write initial rho 
			#}
			
		result$rho_plus[1] <- parms_temp$rho_plus

	result$q_plus_plus  <- vector("numeric", length = length(result$time))		
		result$q_plus_plus[1] <- parms_temp$q_plus_plus
						
			#for(j in 1:length(states)) {
			#result$q_[[1]] <- numeric(length = timesteps/delta+1)
			#result$q_[[1]][1]  <- localdens(x_0, "+", "+")# write initial local densities 
			#}

	result$snapshots <- data.frame(
			ID = parms_temp$ID, 
			assoc = parms_temp$assoc, 
			stock = parms_temp$stock, 
			timestep = 1, 
			g = parms_temp$g, 
			m0 = parms_temp$m0 , 
			b = parms_temp$b, 
			rho_plus = parms_temp$rho_plus, 
			q_plus_plus = parms_temp$q_plus_plus
			)

	result$grids <- list()
	result$grids[[1]] <- initial
	
	result$patches <- list()
	result$patches[[1]] <- patches(initial,"+") 
	
	x_old <- initial
	
for(i in 2:(timesteps/delta+1)) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		if(parms_temp$rho_plus == 0) {flag <- FALSE} else {flag <- TRUE}
		
	 # count local density of occupied fields for each cell: 
		parms_temp$Q_plus <- count(x_old, "+", "+")/4

# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
	if(flag) { # if density is unequal 0, then

		# calculate recolonisation rates of all cells
		recolonisation <- with(parms_temp, (del*rho_plus+(1-del)*Q_plus)*(b-c_*rho_plus)*delta)
		
		# calculate death rates
			# switch for livestock model. 
			rho_comp <- c(0.5, parms_temp$rho_plus)[parms_temp$stock+1]
			
			# for differentiated grazing, i.e. associative protection, 
		if(parms_temp$assoc == TRUE) { 
				# determine vulnerability
				parms_temp$vul <- sum((1-parms_temp$Q_plus)[x_new$cells == "+" & (1-parms_temp$Q_plus) > 0]) /  sum(x_new$cells == "+")#/(width*height)

				death <- with(parms_temp, (m0+g/rho_comp*(1-Q_plus)/vul)*delta)
			
		} else {
			# for undifferentiated grazing

				death <- with(parms_temp, (m0+g/rho_comp)*delta)				
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

		parms_temp$rho_plus <- sum(x_new$cells == "+")/(width*height) 
		result$rho_plus[i] <- parms_temp$rho_plus
		#result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)/delta
		#result$mortality_border[i] <- mean(death[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0], na.rm = TRUE)/delta
		
		parms_temp$rho_plus_plus <- mean(parms_temp$Q_plus[x_new$cells == "+"]/4) 
		result$rho_plus_plus[i] <- parms_temp$rho_plus_plus
		#for(j in 1:length(states)) {
		#	result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
		#	}


		#result$q_[[1]][i]  <- localdens(x_new, states[j], "+")
		

		if(i %in% snapshots ) {
				result$snapshots <- rbind(result$snapshots, data.frame(
										ID = iteration, 
										assoc = parms_temp$assoc, 
										stock = parms_temp$stock, 
										timestep = (i-1)*delta, 
										g = parms_temp$g, 
										m0 = parms_temp$m , 
										b = parms_temp$b, 
										rho_plus =  parms_temp$rho_plus, 
										q_plus_plus = parms_temp$rho_plus_plus
										)
								)

				result$grids[[length(result$grids)+1]] <- x_new 
				result$patches[[length(result$patches)+1]] <- patches(x_new,"+")
						
		}

		# activate to save timeseries at activation of herbivory
		#if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		#if(i %in% seq(1,timesteps, 50)/delta) write.table(c(iteration, i*delta, flag, x_new$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		x_old <- x_new
		
	} # end of simulation.

	result$out <- data.frame(
			ID = iteration,
			starting = parms_temp$starting, 
			assoc = parms_temp$assoc, 
			stock = parms_temp$stock,
			g = parms_temp$g,
			b = parms_temp$b, 
			m0 = parms_temp$m, 			
			#mortality = mean(result$mortality[t_eval], na.rm = TRUE),
			#mortality_border = mean(result$mortality_border[t_eval], na.rm = TRUE),
			rho_plus = mean(result$rho_plus[t_eval], na.rm = TRUE), 
			rho_plus_sd = sd(result$rho_plus[t_eval], na.rm = TRUE), 
			rho_plus_ini = mean(initial$cells == "+"),   #sum(x_0$cells == "+")/(width*height)		
			rho_plus_fin = parms_temp$rho_plus,   #sum(x_0$cells == "+")/(width*height)		
			stable1 = mean(result$rho_plus[t_eval[1:((length(t_eval)-1)/2)] ]  ) - mean(result$rho_plus[t_eval[((length(t_eval)+1)/2):length(t_eval)] ] ) < 0.001,
			stable2 = mean(result$rho_plus[t_eval[1:((length(t_eval)-1)/2)] ]  ) - mean(result$rho_plus[t_eval[((length(t_eval)+1)/2):length(t_eval)] ] ) < 0.01
	)
	
	save(result, file = paste("results/result", iteration, sep = "_"))

} 








	
	##### log-binned patch-size distribution

	if(FALSE) {
	logbins <-  10^seq(log10(1),log10(10000), .1)
	result$logbin <- data.frame(size = logbins)
	
	for(j in 2:length(result$patches)) {
		result$logbin[paste(((c(NA,snapshots)-1)*delta)[j])] <- sapply(1:length(logbins), function(k) length(which(result$patches[[j]] >= logbins[k] & result$patches[[j]] < logbins[k+1])) )
	}
	
	result$logbin$mean <- apply(result$logbin[,2:length(result$patches)], 1, mean)
	result$logbin$sd <- apply(result$logbin[,2:length(result$patches)], 1, sd)
	
	}
	
	
	
	
	#### generating cumulative patch size distribution
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




	#### fitting power-law models to cumulative patch size distribution

result$fit <- list()

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

	# save result to file
	save(result, file = paste("results/result", iteration, sep = "_"))
	
#write.table(out, "output.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
#write.table(result$timeseries, "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)


gc() 
}







