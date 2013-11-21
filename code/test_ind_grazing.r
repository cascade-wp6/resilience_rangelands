
# defining parameter set
parameters = data.frame(
	ID = 1, 
	global = FALSE,
	stock =  TRUE,
	m = .1, 		# intrinsic mortality
	g = .05,		#grazing
	b = .7, 		# beta*eps 
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
prob = c(6/10,2/10,2/10)
color <- c("black","grey80", "white") # define colors for the cell state levels


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

	
# time and resolution of simulation
timesteps = 100
addgrazing = 100
delta = 1/5

t_eval <- ((timesteps/delta)-100/delta):(timesteps/delta)+1
	
snapshots <- c(1, 500/delta+1) #c(addgrazing, timesteps-20:1*50, timesteps)/delta+1)

########################################################################################
########################################################################################
	

highlight <- function(x, colrange = colorRampPalette(c("black", "red2"), space = "rgb"), range = "auto"){

if(is.null(steps)) steps = length(unique(x))
if(range[1] == "auto") {
	min_val = min(x, na.rm = TRUE)
	max_val = max(x, na.rm = TRUE)
	} else {
	min_val = range[1]
	max_val = range[2]
	}


cols <- colrange
steps = length(colrange)
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


greyscale <- colorRampPalette("black","white")

# plotting function for objects of class "landscape".
plot.landscape <- function(x, grid = FALSE, axis = FALSE, cols = "auto", add = FALSE, ani = FALSE, ...) {
  lvls <- levels(x$cells) 
  nlev <- length(lvls)
  if(cols[1] == "auto") cols = color # default color value
  
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

grazingscale <- colorRampPalette(c("#FFFFFF", "#FF0000"), space = "Lab")

points.grazer <- function(g_, cols = "white", add = FALSE, ani = FALSE, ...) {
	
		pos_row <- ceiling((g_$pos)/g_$dim[1])
		pos_col <- g_$pos-g_$dim[1]*(pos_row-1)
		points(pos_row~pos_col, pch = 21, bg = cols, ...)

	
	}
	
plot.grazer <- function(g_, type = "pos", cols = "auto", add = FALSE, ani = FALSE, ...) {
	if(type == "pos") {
		pos_row <- ceiling((g_$pos)/g_$dim[1])
		pos_col <- g_$pos-g_$dim[1]*(pos_row-1)
		points(pos_row~pos_col, pch = 21, bg = col)
	}
	if(type == "dens") {

	#if(cols[1] == "auto") cols =  highlight(g_$cells, grazingscale(16), steps = 16, range = "auto") # default color value
	#if(add) cols = paste(cols, "60", sep = "")
	
#	if(add & cols[1] == "auto")
 cols =  highlight(g_$cells, paste(grazingscale(16), rev(c("F","E","D","C","B","A","9","8","7","6","5","4","3","2","1","0")) , "0", sep = ""),  range = "auto")
  
  if(ani & Sys.info()[['sysname']] == "Windows") adj = -0.5 else adj = 0 #this adjustment constant is added when producing a pixel accurate png or gif file, requires TRUE when the function is used to plot animated figures. 
  
  if(!add) plot(NA,NA, xlim = c(0.5+adj, g_$dim[1]+0.5+adj), ylim = c(g_$dim[2]+0.5+adj, 0+0.5+adj), bty = "n", xaxs = "i", yaxs = "i",xlab = "", ylab = "", xaxt = "n", yaxt = "n", ... ) 

  rect(rep(1:g_$dim[1], times = g_$dim[2])-.5, rep(1:g_$dim[2], each = g_$dim[1])-.5, rep(1:g_$dim[1], times = g_$dim[2])+.5, rep(1:g_$dim[2], each = g_$dim[1])+.5, col = cols, border = NA)
  
	
	}
	
} 

########################################################################################
########################################################################################
	

set.seed(23)
#### first run with grazing
# sampling the initial random grid into a list object
initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = prob ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)

# set position for grazer
n = 1000
steps = 100
delta_t = 0.01

#goat_shack <- matrix(1:prod(x_new$dim), ncol = x_new$dim[1], byrow = TRUE)[(x_new$dim[1]/2-4):(x_new$dim[1]/2+10),(x_new$dim[2]/2-10):(x_new$dim[2]/2+4)]

starting <- sample( which(initial$cells %in% c("0", "-") & count(initial,  c("0", "-")) > 3 ) , n, replace = TRUE) #  & 1:prod(x_new$dim) %in% goat_shack)

g_new <- list()
g_new$dim = initial$dim
g_new$pos = starting
g_new$cells <-  rep(0, length = prod(g_new$dim))
g_new$cells[starting] <- g_new$cells[starting]+delta_t
class(g_new) <- c("list", "grazer")

		
par(mfrow = c(1,1), mar = c(0,0,0,0))

plot(initial, cols = color)
points(g_new)

for(j in 1:steps) {
g_old <- g_new

pos_new_with_border <- sapply(g_old$pos, function(pos) sample(x_to_evaluate[pos]+interact, 1, prob = opencells_logical_with_border[x_to_evaluate[pos]+interact]) )
g_new$pos <- x_with_border[pos_new_with_border]

g_new$cells[g_new$pos] <- g_new$cells[g_new$pos]+delta_t

}


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
	
			result$timeseries <- list()
			result$timeseries$landscape <- list()
			result$timeseries$landscape[[1]] <- initial
			
			
			result$timeseries$grazerdens <- list()
			result$timeseries$grazerdens[[1]] <- g_new
#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	
	parms_temp <- as.list(parameters)
		# copying parms for the simulation 
	#write.table(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, i$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
	#result$timeseries <- as.data.frame(matrix(c( parms_temp$ID, TRUE, FALSE, parms_temp$starting, 1, parms_temp$g, parms_temp$m, parms_temp$b, initial$cells), nrow = 1, dimnames = list("1", header_grids)))

	#parms_temp <- parameters_ini[[c(2,1)[parameters[iteration,4]]]]
	
	#result$patches <- list()
	iteration = 1
for(i in 2:(timesteps/delta+1)) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 
		
		if(parms_temp$rho == 0) {flag <- FALSE} else {flag <- TRUE}
		
		
	if(flag) {
	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		
	# determine vulnerability of a random cell (summed unocupied neighbors over lattice size) 

		parms_temp$dens <- numeric(length = prod(g_new$dim))
		g_with_border <- (g_new$cells)[x_with_border]
			for(k in interact) {
				parms_temp$dens <- parms_temp$dens + g_with_border[x_to_evaluate+k]
			}
		
		parms_temp$vul <- sum(g_new$dens[parms_temp$Q_plus < 1 & parms_temp$Q_plus > 0 ]) # x_old$cells == "+"   g_new$cells
		if(parms_temp$vul == 0) parms_temp$vul = .1
		
		
	
		}
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
	if(flag) { 
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		} else {
		recolonisation <- 0
		}
		
		degradation <- with(parms_temp, (d *delta))
	
	if(flag) {
		if(parms_temp$global == FALSE) {
				death <- with(parms_temp, (m+g*dens/vul)*delta)
		} else {
				death <- with(parms_temp, (m+g/vul)*delta)				
		}
		death[death > 1] <- 1
	}
	
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
				
		#if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		#if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
	if(flag) {		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"}
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
	if(flag) {		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"}
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 run individual grazer dispersion 		
		
starting <- sample( which(x_new$cells %in% c("0", "-") & count(x_new,  c("0", "-")) > 3 ) , n, replace = TRUE) #& 1:prod(x_new$dim) %in% goat_shack

g_new <- list()
g_new$dim = x_new$dim
g_new$pos = starting
g_new$cells <-  rep(0, length = prod(g_new$dim))
g_new$cells[starting] <- g_new$cells[starting]+delta_t
class(g_new) <- c("list", "grazer")

opencells_logical_with_border <- (x_new$cells %in% c("0","-"))[x_with_border]

	
for(j in 1:steps) {
g_old <- g_new

pos_new_with_border <- sapply(g_old$pos, function(pos) sample(x_to_evaluate[pos]+interact, 1, prob = opencells_logical_with_border[x_to_evaluate[pos]+interact]) )
g_new$pos <- x_with_border[pos_new_with_border]

g_new$cells[g_new$pos] <- g_new$cells[g_new$pos]+delta_t

}

# 6 saving state of the new grid		

		result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)/delta
		
		result$mortality_border[i] <- mean(death[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0], na.rm = TRUE)/delta
		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}

		
		for(j in 1:length(states)) {
		result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		}
		
		# activate to save timeseries at activation of herbivory
		#if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		#if(i %in% seq(1,timesteps, 50)/delta) write.table(c(iteration, i*delta, flag, x_new$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		x_old <- x_new
		
		
		result$timeseries$landscape[[i]] <- x_new
		result$timeseries$grazerdens[[i]] <- g_new

	} # end of simulation.

plot(x_new)


plot(x_new, cols = color, add = TRUE)
plot(g_new, type = "dens", add = TRUE)

library("animation")
setwd("C:\\Users\\SCHNEIDER\\Documents\\projects\\CASCADE\\2013 CASCADE grazing\\")
saveGIF( 
	for(i in seq(1, length(result$timeseries$landscape), 1)) {
	  	par(mar = c(0,0,0,0))
			
plot(result$timeseries$landscape[[i]], cols = color, add = FALSE, ani = TRUE)
plot(result$timeseries$grazerdens[[i]], type = "dens", add = TRUE, ani = TRUE)

	}
, movie.name = "ind_based_grazing_23_highb.gif", img.name = "grid", convert = "convert", interval = 0.01/1,
    cmd.fun = system, clean = TRUE, ani.width = width, ani.height = height, outdir = getwd())
	
saveSWF( 
	for(i in seq(1, length(result$timeseries$landscape), 1)) {
	  	par(mar = c(0,0,0,0))
			
plot(result$timeseries$landscape[[i]], cols = color, add = FALSE, ani = TRUE)
plot(result$timeseries$grazerdens[[i]], type = "dens", add = TRUE, ani = TRUE)

	}
, swf.name = "ind_based_grazing_23_highb.swf", img.name = "grid", convert = "convert", interval = 0.01/1,
    cmd.fun = system, clean = TRUE, ani.width = width, ani.height = height, outdir = getwd())
	