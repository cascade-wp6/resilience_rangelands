#################################################
#
#		model: cellular automata helper functions
#
#
#author: Flo
#date: 02.04.2013
#
#
#	Tasks for Flo: 
#		benchmarking tool
#		implement test for stable equilibrium
#		helper functions: randomnumber, count function
# 		plotting functions 
# 		conut patch size and distributions
# 		
#
#################################################


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

greyscale <- colorRampPalette(c( "black", "white"), space = "rgb")

count  <- function(x, neighbor, state = NULL) {
			if(is.null(state)) state = neighbor[1]
			
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
			
			for(k in interact) {
				neighbors <- neighbors + x_logical_with_border[x_state+k]
			}
			# derive a vector which contains the value to count (also "state") in one of the neighbor cells j
			# cumulate the result for each of the four neighboring cells (for(k in interact))

			mean(neighbors/4, na.rm = TRUE)	
			# devide by 4 to get density
			# average over all cells
}


# get patch size and patchsize distribution
patches <- function(x, state) {
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
	
	out <- list()
	out <- sort(patchvec)
	#out <- data.frame(size = unique(out), freq = sapply(unique(out), function(i) length(which(out >= i)) ))
	return(out)
	
	} 

## benchmarking and bottleneck analyses
if(FALSE) {
library(rbenchmark)
benchmark(
 		Qsave = {parms_temp$Q_plus <- rowSums( sapply(interact, 	function(j) (x_old$cells == "+")[x_with_border][x_to_evaluate+j]  ))/4     },
		Qsave_count = {parms_temp$Q_plus <- 	count(x_old, "+")/4     },
		randomnum = {rnum <- runif(width*height)},
		
		rules = {recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, 1-(d *delta))
		death <- with(parms_temp, m*delta)
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		}, 
		applyrules = {
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"
		}, 

# 5 saving state of the new grid		
		savenewrho = {for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}
		},
		savenewq_2 = {
		for(j in 1:length(states)) {
			neighbors <- numeric(length = length(which(x_new$cells == states[j])))
			for(k in interact) {
				neighbors <- neighbors + (x_new$cells == states[j])[x_with_border][x_to_evaluate[x_new$cells == states[j]]+k]
			}
			result$q_[[j]][i]  <- mean(neighbors/4)	
		}
		},
		
		savenewq_3 = {
		for(j in 1:length(states)) {
			result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		}
		
		}, replications = 5000)

		}
		
		
