#################################################
#
#		model: iterator for cellular automata
#
#
#author: Flo
#date: 02.04.2013
#
#
#
#################################################



rm(list=ls())
source("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\automatafunctions.r")
setwd()
library(foreach)
library(doSNOW)

cl <- makeCluster(rep("localhost", times = 11), type = "SOCK")
registerDoSNOW(cl)

# specify lattice
width = 50
height = 50


# time and resolution of simulation
timesteps = 10
addgrazing = 1
delta = 1/10


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


# defining parameter set
parameters = list(
	m = 0.1, 		# intrinsic mortality
	b = 0.5, 		# beta*eps 
	d = 0.2,		# degradation
	c_ = 0.0, 		# beta*g  
	del = 0.0,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.0, 	# regeneration rate
	f = 0.8, 		# local fascilitation
	g = 0.05			# grazing per neighboring unvegetated cell
)
	
# initial cell states
states = c("+","0","-")
color <- c("black","grey80", "white") # define colors for the cell state levels
	
env <- 	seq(0.0,1,length = 50)
grazing <- seq(0, 0.25, length = 3)



lgrazing <- length(grazing)
lenv <- length(env)
env <- rep(env, times = lgrazing)
grazing <- rep(grazing, each = lenv)

foreach(iteration = 1:(lenv*lgrazing)) %dopar% {

parameters$b <- env[iteration]
parameters$g <- grazing[iteration]

flag = "first"

while(flag %in% c("first", "checkbistability")) {

set.seed(iteration)
if(flag == "first") prob = c(.1/10,.9/10,9/10)
if(flag == "checkbistability") prob = c(8/10,.9/10,1.1/10)

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
			
			#save all level global densities: rho_*
			result$rho  <- list()			
			for(i in 1:length(states)) {
				result$rho[[i]] <- numeric(length = timesteps/delta)
				result$rho[[i]][1]  <- sum(initial$cells == states[i])/(width*height) # write initial rho 
			}
			
			result$q_  <- list()			
			for(j in 1:length(states)) {
			result$q_[[j]] <- numeric(length = timesteps/delta)
			result$q_[[j]][1]  <- mean(rowSums( sapply(interact, 	function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == states[j]]+k]) )) /4# write initial local densities 
			}
			
			result$timeseries <- list()
			result$timeseries$initial <- initial
			result$timeseries$stable_without_grazing <- initial
			result$timeseries$final_with_grazing <- initial
			
#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	parms_temp <- parameters # copying parms for the simulation and multiplying
	
for(i in 1:(timesteps/delta)+1) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 
		
	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		# This is a vectorised evaluation: first transform the original grid x_old, into a grid with border. Here check for all cells of the original grid (x_to_evaluate) the content of all four neighbouring cells, whose relative position is defined in interact. The number of cells with content "+" are count and divided by 4.
	# count local density of degraded fields for each cell
		parms_temp$Q_unveg <- (count(x_old, "0")+count(x_old, "-"))/4 
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		if(i < addgrazing/delta) {death <- with(parms_temp, m*delta)} else {death <- with(parms_temp, (m+g*Q_unveg)*delta)}
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > .95 )) warning(paste("a set probability is exceeding 1 in time step", i, "! decrease delta!!!")) 
				
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}

		
		for(j in 1:length(states)) {
		result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		}
		
		# activate to save each single timeseries step
		if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		
		# activate for plotting during simulation (very slow !!!)
		#if(i %in% seq(1,(timesteps/delta)+1, 10)) plot(x_new, col = color, add = TRUE)

		x_old <- x_new
		
		gc()  #garbage collection
	} # end of simulation.

	result$timeseries$final_with_grazing <- x_new
gc()
	t_eval <- ((timesteps/delta)-50/delta):(timesteps/delta)+1
	
	out <- c(grazing = parameters$g, 
			environmt = parameters$b,  
			rho_plus = mean(result$rho[[1]][t_eval]), 
			rho_plus_ungrazed = mean(result$rho[[1]][(t_eval-1)-(timesteps-addgrazing)/delta]), 
			rho_minus = mean(result$rho[[2]][t_eval]),
			q_plus = mean(result$q_[[2]][t_eval]))
			
out <- list(grazing = parameters$g, 
			environmt = parameters$b,  
			rho_plus = mean(result$rho[[1]][t_eval]), 
			rho_plus_ungrazed = mean(result$rho[[1]][(t_eval-1)-(timesteps-addgrazing)/delta]), 
			rho_minus = mean(result$rho[[2]][t_eval]),
			q_plus = mean(result$q_[[2]][t_eval]),
			t_ini = result$timeseries$initial,
			t_stb = result$timeseries$stable_without_grazing,
			t_fin = result$timeseries$final_with_grazing)
			
if(flag == "checkbistability") {flag = "ok"}
if(flag == "first" && out$rho_plus < 10e-4) {flag = "crashed"} else {flag = "ok"}


#}

out
}
} -> saveall


output <- as.data.frame(saveall[[1]][1:4])
for(i in 2:length(saveall) ) {
output <- rbind(output, saveall[[i]][1:4])
}


par(mfrow = c(3,1))
with(data = output[output$grazing == 0,], plot(environmt, rho_plus, ylim = c(0,1),type = "b"))
with(data = output[output$grazing == 0.125,], plot(environmt, rho_plus, ylim = c(0,1),type = "b"))
with(data = output[output$grazing == 0.25,], plot(environmt, rho_plus, ylim = c(0,1),type = "b"))


par(mfrow = c(1,2), mar = c(1,1,1,1))
plot(saveall[[32]]$t_ini, col = color)
plot(saveall[[32]]$t_fin, col = color)



stopCluster(cl)
