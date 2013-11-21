setwd("C:\\Users\\SCHNEIDER\\Documents\\CASCADE\\2013_grazingmodel\\corrected_grazing_vul")

env <- 	seq(0.2,1,.04)#seq(0.2,1,.025)
graz <- seq(0,.1, 0.005)#seq(0,.9, 0.075)
mort <- 0.1 #seq(0.1,.3, 0.1)
init <- c("low", "high")
global <- c(TRUE, FALSE)

lgraz <- length(graz)
lenv <- length(env)
lmort <- length(mort)
linit <- 2
lglobal <- 2

# defining parameter set
parameters = data.frame(
	ID = 1:(lgraz*lmort*lglobal*linit*lenv), 
	global = rep(global, each = lgraz*lenv*lmort*linit),
	starting = rep(rep(init, each = lgraz*lenv*lmort), times = lglobal),
	g = rep(rep(graz, each = lenv), times = lmort*lglobal*linit),		#grazing
	m = rep(rep(mort, each = lgraz*lenv), times = lglobal*linit), 		# intrinsic mortality
	b = rep(env, times = lgraz*lmort*lglobal*linit), 		# beta*eps 
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
  if(cols[1] == "auto") cols = greyscale(nlev) # default color value
  
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

	
output <- read.csv("output.csv")

	

plotgrid <- function(x , ...) { 
load(paste("results\\result", x, sep = "_"))

grd <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = as.factor(states[as.integer(result$timeseries[23,-(1:7)])])  #second element contains a random row-wise, factorial vector to fill the grid 
	)
class(grd) <- c("list","landscape") # set class of object (required for plotting)

plot(grd, cols = color[match(levels(grd$cells), states)], ...)

}

	load(paste("results\\result", 468, sep = "_"))

grd <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = as.factor(states[as.integer(result$timeseries[23,-(1:7)])])  #second element contains a random row-wise, factorial vector to fill the grid 
	)
class(grd) <- c("list","landscape") # set class of object (required for plotting)

plot(grd, cols = color[match(levels(grd$cells), states)], ...)

	
patches(grd, "+")


plotgrid(468)
box()




load(paste("results\\result", 22, sep = "_"))


filenames <- list.files("results\\")


library(foreach)
library(doSNOW)

workerlist <- rep("localhost", times = 11)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)

foreach(i = filenames) %dopar% {
load(paste("results\\", i, sep = ""))

	#bins <- sort(unique(unlist(result$patches[-1])))
	bins <- 10^seq(log10(1),log10(10000), .1)
	#bins <- seq(1,10001, 10)

	if( !any(is.na(unlist(result$patch[3:22])))) {
	result$cumpatch2 <- data.frame(size = bins)
	
	for(j in 2:length(result$patches)) {
		result$cumpatch2[paste(((snapshots-1)*delta)[-1][j])] <- sapply(bins, function(k) length(which(result$patches[[j]] >= k)) )
	}
	
	result$cumpatch2$mean <- apply(result$cumpatch2[,2:length(result$patches)], 1, mean)
	result$cumpatch2$sd <- apply(result$cumpatch2[,2:length(result$patches)], 1, sd)
	
	result$cumpatch2$mean_uni <- result$cumpatch2$mean/max(result$cumpatch$mean)
	result$cumpatch2$sd_uni <- result$cumpatch2$sd/max(result$cumpatch$mean)
	
	} else {
	result$cumpatch2 <- NA
	}


save(result, file = (paste("results\\", i, sep = "")) )

#file.remove(paste("results\\", i, sep = ""))
}


stopCluster(cl)



