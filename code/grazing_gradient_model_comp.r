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
	

highlight <- function(x, colrange = c("black", "red2"), steps = NULL, range = "auto"){

if(is.null(steps)) steps = length(unique(x))
if(range == "auto") {
	min_val = min(x, na.rm = TRUE)
	max_val = max(x, na.rm = TRUE)
	} else {
	min_val = range[1]
	max_val = range[2]
	}

colorscale <- colorRampPalette(colrange, space = "rgb")
cols <- colorscale( steps)

cols[as.integer((x-min_val)/(max_val-min_val)*(steps-1)+1)]

}

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
	old <- rep(1, times = prod(x$dim)) 
	
	while(!identical(old[pattern], map[pattern])) {
		old <- map
		count = as.integer(1)
		for(i in which(pattern)) {
			if(all(is.na(map[x_with_border][x_to_evaluate[i]+interact])) ) {map[i] <- count}
			count <- count +1
			if(any(!is.na(map[x_with_border][x_to_evaluate[i]+interact])) ) map[i] <- min(map[x_with_border][x_to_evaluate[i]+interact], na.rm = TRUE)
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
first_ID = 1
env <- 	seq(0.2,1,.025) #seq(0.2,1,.025)
graz <- seq(.0, .1, 0.025) #.025 #
mort <- 0.05 #c(0.05, 0.1) #seq(0.05,.15, 0.05)
init <- c(0.25, 0.9)  #seq(0.25, 1, .05)
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

t_eval <- ((timesteps/delta)-100/delta):(timesteps/delta)+1
	
	

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

out_header <- data.frame(
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
			rho_plus_ini = NA,
			q_plus = NA,
			rho_zero = NA,
			rho_minus = NA
			)

write.table(out_header[-1,], "output.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)

header_grids <- c("ID", "global", "stock", "starting", "timestep", "g", "m0" , "b", "rho_plus", "q_plus", paste("X", 1:(width*height), sep =""))

dir.create("results")

################ starting parallel backend
	
library(foreach)
library(doSNOW)

#workerlist <- rep(list(ubuWorker), times = 2)
workerlist <- rep("localhost", times = 23)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)


################ starting foreach loop
		
foreach(iteration = parameters$ID) %dopar% {

#iteration = 22
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
		b = parameters[iteration,4] , 		# beta*eps 
		d = 0.1,		# degradation
		c_ = 0.2, 		# beta*g  
		del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
		r = 0.01, 	# regeneration rate
		f = 0.9, 		# local fascilitation
		g = 0		#grazing
		)
	
for(i in 1:(1000/delta+1)) {    #calculation loop
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
				result$rho[[i]][1]  <- sum(initial$cells == states[i])/(width*height) # write initial rho 
			}
			
			result$q_  <- list()			
			for(j in 1:length(states)) {
			result$q_[[j]] <- numeric(length = timesteps/delta+1)
			result$q_[[j]][1]  <- mean(rowSums( sapply(interact, function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == "+"]+k]) )) /4# write initial local densities 
			}

#  initialising first iteration 
	x_old <- x_0    #old landscape for first timestep   
	
	parms_temp <- as.list(parameters[iteration,])
		# copying parms for the simulation 
	#write.table(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, i$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
	result$timeseries <- as.data.frame(matrix(c( parms_temp$ID, TRUE, FALSE, parms_temp$starting, 1, parms_temp$g, parms_temp$m, parms_temp$b, result$rho[[1]][1], result$q_[[1]][1], initial$cells), nrow = 1, dimnames = list("1", header_grids)))

	result$patches <- list()
	result$patches[[1]] <- patches(x_0,"+") 
	
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

		
		for(j in 1:length(states)) {
		result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		}
		
		if(i %in% snapshots) {
				result$timeseries <- rbind(result$timeseries, 
								as.data.frame(matrix(c( iteration, parms_temp$global, parms_temp$stock, parms_temp$starting, (i-1)*delta, parms_temp$g, parms_temp$m, parms_temp$b,  result$rho[[1]][i], result$q_[[1]][i],  x_new$cells), nrow = 1, dimnames = list("1", header_grids)))
				)
				result$patches[[length(result$patches)+1]] <- patches(x_new,"+") 
		}

		# activate to save timeseries at activation of herbivory
		#if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		#if(i %in% seq(1,timesteps, 50)/delta) write.table(c(iteration, i*delta, flag, x_new$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		x_old <- x_new
		
	} # end of simulation.
	
	
	
	#bins <- sort(unique(unlist(result$patches[-1])))
	bins <- 10^seq(log10(1),log10(10000), .1)


	if(length(bins) > 0) {
	result$cumpatch <- data.frame(size = bins)
	
	for(i in 2:length(result$patches)-1) {
		result$cumpatch[paste(((snapshots-1)*delta)[i])] <- sapply(bins, function(j) length(which(result$patches[[i]] >= j)) )
	}
	
	result$cumpatch$mean <- apply(result$cumpatch[,2:length(result$patches)], 1, mean)
	result$cumpatch$sd <- apply(result$cumpatch[,2:length(result$patches)], 1, sd)
	} else {
	result$cumpatch <- NA
	}

	# save result to file
	save(result, file = paste("results/result", iteration, sep = "_"))
	
	out <- data.frame(
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
			rho_plus_ini = mean(x_0$cells == "+"),   #sum(x_0$cells == "+")/(width*height) 
			q_plus = mean(result$q_[[1]][t_eval], na.rm = TRUE),
			rho_zero = mean(result$rho[[2]][t_eval], na.rm = TRUE),
			rho_minus = mean(result$rho[[3]][t_eval], na.rm = TRUE)
			)

write.table(out, "output.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)


#write.table(result$timeseries, "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)


gc() 
}


stopCluster(cl)



filenames <- list.files("results/")

output <- read.csv("output.csv")


		
		