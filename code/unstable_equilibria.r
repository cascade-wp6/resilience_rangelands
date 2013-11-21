

###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
# unstable equilibria
# 



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

# specify lattice
width = 100
height = 100

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
#output_fits <- read.csv("output_fits.csv")

output <- output[ output$g != 1 & output$stable == TRUE, ] # & output$ID > 27541

# initial cell states
states = c("+","0","-")
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
timesteps = 2000
delta = 1/5


library(foreach)
library(doSNOW)

#workerlist <- rep(list(ubuWorker), times = 2)
workerlist <- rep("localhost", times = 15)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)



		out <- data.frame(	global =  NA,
			stock =  NA,
			m = 0.05, 		# intrinsic mortality
			g = NA,		#grazing
			b = NA, 		# beta*eps 
			d = 0.1,		# degradation
			c_ = 0.2, 		# beta*g  
			del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
			r = 0.01, 	# regeneration rate
			f = 0.9,		# local fascilitation
			rho_plus = NA,
			vul = NA,
			eq = NA
			)[-1,]
			
write.table(out, "unstable_eq.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)

for(g in seq(0,.1,0.025)) {
	for(l in 1:4) {	

		foreach(b = seq(0.2,1, 0.005), .combine = "rbind") %dopar% {
		#foreach(b = seq(0.2,1, 0.01), .combine = "rbind") %do% {
		#for(b in seq(0.2,1, 0.01)) {

		take = output[which(as.character(output$g) == as.character(g) & output$b < b+0.01 & output$b > b-0.01 & output$stock == c(FALSE, FALSE, TRUE, TRUE)[l] & output$globalgrazing == c(TRUE, FALSE, TRUE, FALSE)[l] & output$stable == TRUE),]
		
		
		high = max(take$rho_plus)
		low = 0
		eq = mean(c(high, low))

		for(i in 1:15) {
			
		initial <- list(  
			dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
			cells = sample(factor(1:length(states)), width*height, replace = T, prob = c(eq, (1-eq)/2, (1-eq)/2) ) #second element contains a random row-wise, factorial vector to fill the grid 
		)
		levels(initial$cells) <- states  #assign cell states 
		class(initial) <- c("list","landscape") # set class of object (required for plotting)

		#plot(initial, ani = TRUE)
	
		x_old <- initial
		
		parms_temp <- list( 
			global =  c(TRUE, FALSE, TRUE, FALSE)[l],
			stock =  c(FALSE, FALSE, TRUE, TRUE)[l],
			m = 0.05, 		# intrinsic mortality
			g = g,		#grazing
			b = b, 		# beta*eps 
			d = 0.1,		# degradation
			c_ = 0.2, 		# beta*g  
			del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
			r = 0.01, 	# regeneration rate
			f = 0.9,		# local fascilitation
			rho = sum(x_old$cells == "+")/(width*height)  ,
			vul = NA,
			equlbr = NA
		)
		
		
		for(j in 2:101) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		
		if(parms_temp$rho == 0) {flag <- FALSE} else {flag <- TRUE}
		
		
	 # count local density of occupied fields for each cell: 
		Q_plus <- count(x_old, "+")/4   

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
				parms_temp$vul <- sum((1-Q_plus)[x_new$cells == "+" & (1-Q_plus) > 0]) /  sum(x_new$cells == "+")#/(width*height)

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
		#if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		#if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		# apply rules 
	if(flag) {	
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		}
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"

		#rho_t <- cbind(rho_t, parms_temp$rho)
		x_old <- x_new
		parms_temp$rho <- sum(x_new$cells == "+")/(width*height)   # get old rho plus for this timestep 

}
	
	if(parms_temp$rho >= eq ) {
			high = eq
			eq = mean(c(eq, low))
	} else {
			low = eq
			eq = mean(c(eq, high))
	}
		
		
		
		}
		parms_temp$equlbr <- eq
		
		#out <- rbind(out, as.data.frame(parms_temp))
		return(as.data.frame(parms_temp))
		} -> out_temp
		
		#take = output[which(output$g == g & output$stock == c(FALSE, FALSE, TRUE, TRUE)[l] & output$globalgrazing == c(TRUE, FALSE, TRUE, FALSE)[l] & output$stable == TRUE),]
			
		#plot(out_temp$b, out_temp$eq, ylim = c(0,1), pch = 20, col = "gray80")
		#points(take$b, take$rho_plus , pch = 20, col = "gray20")
	
		write.table(out_temp, "unstable_eq.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = TRUE)

		#out <- rbind(out, out_temp)
		
	gc() 

	}
}


stopCluster(cl)
#write.table(out, "unstable_eq.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)

#save(out, file = "unstable_eq")
