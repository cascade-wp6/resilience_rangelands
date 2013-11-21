################################################
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
#source("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\automatafunctions.r")
#setwd("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\thread herbivory\\")

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

library(foreach)
library(doSNOW)

winWorker <- list(host = "localhost", 
				rscript = "C:\\R\\R-2.15.3\\bin\\x64\\Rscript.exe", 
				snowlib = "C:\\R\\R-2.15.3\\library"
	)
winWorker <- list(host = "localhost", 
				rscript = "C:/R/R-2.15.3/bin/x64/Rscript.exe", 
				snowlib = "C:/R/R-2.15.3/library"
	)
ubuWorker <- list(host = "kefi118", user = "schneider",
				rscript = "/usr/lib/R/bin/Rscript", #R/bin/
				snowlib = "/usr/lib/R/site-library/"
	)
#	, rep(list(ubuWorker), times = 23)
workerlist <- rep("localhost", times = 23)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)

# specify lattice
width = 50  
height = 50


# time and resolution of simulation
timesteps = 2000
addgrazing = 1000
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
# Effect of mortality (m). d = 0.1, r = 0.01, c = 0.2, d = 0.1, f = 0.9. (Kefi et al 2007 TPB Fig 4c)
parameters = list(
	m = 0.1, 		# intrinsic mortality
	b = 0.5, 		# beta*eps 
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.01, 	# regeneration rate
	f = 0.9, 		# local fascilitation
	g = 0.025		#grazing
)
	
# initial cell states
states = c("+","0","-")
color <- c("black","grey80", "white") # define colors for the cell state levels
	
env <- 	seq(0.2,1,.1)
graz <- seq(0,.9, 0.05)
mort <- seq(0.1,.5, 0.1)

lgraz <- length(graz)
lenv <- length(env)
lmort <- length(mort)

env <-	rep(rep(env, times = lgraz), times = lmort)
graz <-	 rep(rep(graz, each = lenv), times = lmort)
mort <-  rep(mort, each = lgraz*lenv)

replicates = 2
init <- rep(c("low", "high"), each = length(mort) )
env <- rep(env, times = replicates)
graz <- rep(graz, times = replicates)
mort <- rep(mort, times = replicates)


str(data.frame(init, env, graz, mort))


header_output <- data.frame(
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
			rho_minus = NA)
			

write.table(header_output[-1,], "output.csv", row.names = FALSE, sep = ",")
header_grids <- c("iteration", "i", "g", "m" , "b", paste("X", 1:(width*height), sep =""))

write.table(as.data.frame(matrix(header_grids, nrow = 1, dimnames = list(NA, header_grids)))[-1,], "grids.csv", row.names = FALSE, sep = ",", append = TRUE)

#snapshots <- c(0)
#snapshots <- c(1, c(addgrazing, timesteps)/delta+1)
snapshots <- c(1, c(addgrazing, timesteps-4:1*50, timesteps)/delta+1)

foreach(iteration = 1:length(env)) %dopar% {

parameters$b <- env[iteration]
parameters$g <- graz[iteration]
parameters$m <- mort[iteration]

flag = "high"

while(flag %in% c("low", "high")) {

set.seed(iteration)
if(flag == "low") prob = c(.1/10,.9/10,9/10)
if(flag == "high") prob = c(9/10,.9/10,0.1/10)

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
			
			result$mortality <- numeric(length = timesteps/delta)
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
			
			#result$timeseries <- list()
			#result$timeseries$initial <- initial
			
#write.table(matrix(c(iteration, flag, initial$cells), nrow = 1), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	parms_temp <- parameters # copying parms for the simulation and multiplying
	#write.table(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, i$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
	if(1 %in% snapshots) write.table(as.data.frame(matrix(c(iteration, 1, parms_temp$g, parms_temp$m, parms_temp$b, initial$cells), nrow = 1, dimnames = list(NA, header_grids))), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

	
for(i in 1:(timesteps/delta)+1) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 

	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		# This is a vectorised evaluation: first transform the original grid x_old, into a grid with border. Here check for all cells of the original grid (x_to_evaluate) the content of all four neighbouring cells, whose relative position is defined in interact. The number of cells with content "+" are count and divided by 4.
	# count local density of degraded fields for each cell
		parms_temp$Q_unveg <- ( count(x_old, c("0", "-")))/4 
		parms_temp$vul <- sum(parms_temp$Q_unveg[x_old$cells == "+"]*4)/(width*height)
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		if(i < addgrazing/delta) {death <- with(parms_temp, m*delta)} else {death <- with(parms_temp, (m+g/vul*Q_unveg)*delta)}
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
		
		if(i %in% snapshots) write.table(as.data.frame(matrix(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, x_new$cells), nrow = 1, dimnames = list(NA, header_grids))), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		# activate to save timeseries at activation of herbivory
		#if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		#if(i %in% seq(1,timesteps, 50)/delta) write.table(c(iteration, i*delta, flag, x_new$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		x_old <- x_new
		
		gc()  #garbage collection
	} # end of simulation.

gc()
	t_eval <- ((timesteps/delta)-100/delta):(timesteps/delta)+1
	
	out <- data.frame(
			ID = iteration,
			starting = flag, 
			grazing = parameters$g, 
			environmt = parameters$b, 
			mort_int = parameters$m, 			
			mortality = mean(result$mortality[t_eval]),
			mortality_border = mean(result$mortality_border[t_eval], na.rm = TRUE),
			rho_plus = mean(result$rho[[1]][t_eval]), 
			rho_plus_ungrazed = mean(result$rho[[1]][(t_eval-1)-(timesteps-addgrazing)/delta]), 
			q_plus = mean(result$q_[[1]][t_eval]),
			q_plus_ungrazed = mean(result$q_[[1]][(t_eval-1)-(timesteps-addgrazing)/delta]),
			rho_zero = mean(result$rho[[2]][t_eval]),
			rho_minus = mean(result$rho[[3]][t_eval])
			)

write.table(out, "output.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

if(flag == "high") {flag = "ok"}
if(flag == "low" && out$rho_plus < 10e-3) {flag = "high"} else {flag = "ok"}


#}

}


} 


if(FALSE) {
output <- read.csv("grazing_gradient\\output.csv", header = TRUE)   

highlight <- function(x, colrange = c("black", "red2"), steps = NULL){

if(is.null(steps)) steps = length(unique(x))
 
colorscale <- colorRampPalette(colrange, space = "rgb")
cols <- colorscale( steps)

cols[as.integer((x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))*steps+1)]

}



grazingscale <- colorRampPalette(c("black", "red2"), space = "rgb")
cols <- grazingscale(50)

pdf("grazing_baseline.pdf", height = 8, width = 4, paper = "special")
par(mfrow = c(4,1), mar = c(4,4,1,1))
with(data = output, plot(environmt, rho_plus, xlim = c(0,1), ylim = c(0,1), pch = 20, col = highlight(grazing) ))
with(data = output, plot(environmt, q_plus, xlim = c(0,1), ylim = c(0,1), pch = 20, col = highlight(grazing) ))
with(data = output, plot(environmt, mortality, xlim = c(0,1), ylim = c(0,1),  pch = 20, col = highlight(grazing) ))
with(data = output, plot(environmt, mortality_border, xlim = c(0,1), ylim = c(0,1), pch = 20, col = highlight(grazing) ))

dev.off()

with(data = output, plot(grazing, mortality, ylim = c(0,.4), xlim = c(0,.6), pch = 20, col = highlight(environmt, c("black", "green3")) ))
with(data = output, points(grazing, mortality_border, pch = 20, col = highlight(environmt, c("black", "red3")) ))

pdf("grazing_gradient2.pdf", height = 3, width = 3, paper = "special")
par(mfrow = c(1,1), mar = c(4,4,2,2))
for(i in sort(unique(output$environmt))) with(data = output[output$environmt == i,], plot(rho_plus, mortality,  main = paste("environmt =",i), pch = 20 , col =   highlight(grazing), xlim = c(0,1), ylim = c(0,.9), cex = 3))
dev.off()

with(data = output, plot(rho_plus, mortality, pch = 20 , col =  highlight(environmt, c("black", "beige","blue"))))
with(data = output, plot(rho_plus, mortality, pch = 20 , col =  highlight(grazing)))



with(data = output, plot(mortality, rho_plus_ungrazed, pch = 20 , col =  highlight(environmt, c("black", "beige","blue"))))

with(data = output, plot(rho_plus_ungrazed,rho_plus,  pch = 20 , col =  highlight(environmt, c("black", "beige","blue"))))

with(data = output, plot(rho_plus_ungrazed,rho_plus,  pch = 20 , col =  highlight(mortality_border,  c("black", "red3"))))
par(mfrow = c(3,3))
for(i in sort(unique(output$environmt))) with(data = output[output$environmt == i,], plot(mort_int,grazing,  pch = 20 , col =  highlight(rho_plus,  c("grey95", "darkgreen")), xlim = c(0,1), ylim = c(0,.9), cex = 3))

pdf("grazing_mortality1.pdf", height = 2.5, width = 6, paper = "special")
par(mfrow = c(1,3))
for(i in c(0.1, 0.2, 0.3)) with(data = output[output$mort_int == i,], plot(environmt,grazing,  pch = 20 , main = paste("mortality =",i), col =  highlight(rho_plus,  c("grey95", "darkgreen")), xlim = c(0,1), ylim = c(0,.9), cex = 3))
dev.off()

pdf("grazing_gradient.pdf", height = 3, width = 3, paper = "special")
par(mfrow = c(1,1), mar = c(4,4,2,2))
for(i in sort(unique(output$grazing))) with(data = output[output$grazing == i,], plot(environmt,mort_int,  main = paste("grazing =",i), pch = 20 , col =  highlight(rho_plus,  c("grey95", "darkgreen")), xlim = c(0,1), ylim = c(0,.9), cex = 3))
dev.off()

pdf("grazing_mortality2.pdf", height = 6, width = 3, paper = "special")
par(mfrow = c(3,1))
for(j in sort(unique(output$grazing))) {
for(i in c(0.1, 0.2, 0.3)) with(data = output[output$mort_int == i & output$grazing == j,], plot(environmt,rho_plus,  pch = 20 , main = paste("mortality =",i), xlim = c(0,1), ylim = c(0,.9), cex = 3))


}
dev.off()

}

stopCluster(cl)
