library(rbenchmark)



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


# specify lattice
width = 100
height = 100

# initial cell states
states = c("+","0","-")
prob = c(9/10,.9/10,0.1/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
timesteps = 1500
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

snapshots <-  (timesteps-10:0*50)/delta+1


benchmark(

updates = {

  #calculation loop
for(i in 2:501) {
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

		x_old <- x_new
		}
		
},

savresults = {
	for(i in 2:501) {
# 5 saving state of the new grid		

		result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)/delta
		
		result$mortality_border[i] <- mean(death[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0], na.rm = TRUE)/delta
		
		for(j in 1:(length(states)-1)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}

		
		}
		},
		
savresults2a = {
	for(i in 2:501) {		
		result$q_[[1]][i]  <- localdens(x_new, "+", "+")
		
		}
		},
		
			
savesnapshots1 = {
		 
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
								
		},

savesnapshots2 = {
				result$timeseries[[length(result$timeseries)+1]] <- x_new 
				},
savesnapshots3 = {
				result$patches[[length(result$patches)+1]] <- patches(x_new,"+")
					}	



, replications = 3)



	