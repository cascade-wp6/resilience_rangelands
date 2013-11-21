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
source("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\automatafunctions.r")
source("E:\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\automatafunctions.r")
setwd("E:\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\thread herbivory\\")
setwd("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\thread herbivory\\")

# specify lattice
width = 50
height = 50

# initial cell states
states = c("+","0","-")
prob = c(1/10,1/10,8/10)

# time and resolution of simulation
timesteps = 200
addgrazing = 200
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

	

# sampling the initial random grid into a list object
initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = prob ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)

color <- c("black","grey80", "white") # define colors for the cell state levels

#plot the initial state
par(mar = c(0,0,0,0), mfrow = c(1,1))
plot(initial, cols = color)

# defining parameter set
parameters = list(
	m = 0.1, 		# intrinsic mortality
	b = 0.9, 		# beta*eps 
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.01, 	# regeneration rate
	f = 0.9, 		# local fascilitation
	g = .1		#grazing
)

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
			result$q_[[j]][1]  <- mean(rowSums( sapply(interact, 	function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == states[j]]+k]) )) /4# write initial rho 
			}
			
			result$timeseries <- list() # create a subordinate list and
			result$timeseries[[1]] <- initial # write 'x' as the first entry of the list
			for(i in 1:(timesteps/delta)+1) result$timeseries[[i]] <- initial	# allocate memory for each timeseries object
			

#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	parms_temp <- parameters # copying parms for the simulation and multiplying
	
	
	
for(i in 1:(timesteps/delta)+1) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 
		#parms_temp$vul <- sum(count(x_old, c("0", "-"), "+" )/4)/(width*height)
	
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
		
		if(any(c(recolonisation+degradation, death, regeneration) > .95 )) stop(paste("a set probability is exceeding 1 in time step", i, "! decrease delta!!!")) 
				
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) stop(paste("a set probability falls below 0 in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid	

		result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)
		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}

		
		for(j in 1:length(states)) {
		result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		}
		
		# activate to save each single timeseries step
		result$timeseries[[i]] <- x_new  #the whole grid is saved to timeseries
		
		# activate for plotting during simulation (very slow !!!)
		#if(i %in% seq(1,(timesteps/delta)+1, 10)) plot(x_new, col = color, add = TRUE)

		x_old <- x_new 
		gc()  #garbage collection
	} # end of simulation.


#str(result)

	
highlight <- function(x, colrange = c("black", "red2"), steps = NULL){

if(is.null(steps)) steps = length(unique(x))
 
colorscale <- colorRampPalette(colrange, space = "rgb")
cols <- colorscale( steps)

cols[as.integer((x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))*steps+1)]

}

# the rest is graphical output

# FIGURE 1 	-- final states of the grid

#pdf("grazing_corrected.pdf", width = 12, height = 6, paper = "special" )
windows(10,5)
# FIGURE 2 	-- global densities and local densities around cells in state "+" 
layout(matrix(c(1,1,2,3), ncol = 2))
par(mar = c(2,2,2,2))
plot(x_new, grid = FALSE, cols = color)
box()

par(mar = c(2.5,4,2,1)+.1, xaxt = "s", yaxt ="s")
plot(NA,NA, ylim = c(0, 1), xlim = range(result$time), type ="l", xlab = "time", ylab = expression(rho ['*']), xaxs ="i", yaxs = "i" )
lines(result$time, result$rho[[1]])
polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]]+result$rho[[2]]+result$rho[[3]], 0), col = "white")
polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]]+result$rho[[2]], 0), col = "grey50")
polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]], 0), col = "black")
	
par(mar = c(4.5,4,0,1)+.1, xaxt = "s", yaxt ="s")
plot(NA,NA, ylim = c(0, 1), xlim = range(result$time), type ="l", xlab = "time", ylab = "q*", xaxs ="i", yaxs = "i" )
polygon(c(0,result$time,max(result$time)), c(0, result$mortality, 0), col = "black")
par(mfrow = c(1,1))
	
#dev.off()
	
plot(result$rho[[1]],result$mortality/delta, col = highlight(result$mortality) , xlim = c(0,1), pch = 20 )
par(mfrow = c(1,1))


par(mar = c(4.5,4,0,1)+.1, xaxt = "s", yaxt ="s")
plot(NA,NA, ylim = c(0, 1), xlim = range(result$time), type ="l", xlab = "time", ylab = "q*", xaxs ="i", yaxs = "i" )
lines(result$time, result$q_[[1]])
polygon(c(0,result$time,max(result$time)), c(0, result$q_[[1]]+result$q_[[2]]+result$q_[[3]], 0), col = "white")
polygon(c(0,result$time,max(result$time)), c(0, result$q_[[1]]+result$q_[[2]], 0), col = "grey50")
polygon(c(0,result$time,max(result$time)), c(0, result$q_[[1]], 0), col = "black")
par(mfrow = c(1,1))
dev.off()


par(mar = c(2.5,4,2,1)+.1, xaxt = "s", yaxt ="s")
plot(result$time, result$mortality/delta, type ="l", lwd = 2)


grazingscale <- colorRampPalette(c("black", "red"), space = "rgb")
grazingcol <- c(grazingscale(12), color[2:3] )


patch <- c(patches(result$timeseries[[10001]], "+"),patches(result$timeseries[[9501]], "+"),patches(result$timeseries[[9001]], "+"), patches(result$timeseries[[8501]], "+"), patches(result$timeseries[[8001]], "+"), patches(result$timeseries[[7501]], "+"))



pdf("herbivory_high_patchsizedist.pdf", width = 12, height = 6, paper = "special" )

layout(matrix(c(1,1,2,3), ncol = 2))
par(mar = c(2,2,2,2))
		plot(x, grid = FALSE, cols = grazingcol[as.numeric(levels(x$cells))], ani = FALSE)
		box()
par(mar = c(2.5,4,2,1)+.1, xaxt = "s", yaxt ="s")

plot(data.frame(size = unique(patch), freq = sapply(unique(patch), function(i) length(which(patch >= i)) )), log = "xy", pch = 20)

par(mar = c(4.5,4,0,1)+.1, xaxt = "s", yaxt ="s")

bins <- seq(1,max(patch),5)
plot(data.frame(size = bins, freq = sapply(bins, function(i) length(which(patch >= i & patch < i+5)) )), log = "xy", pch = 20)
dev.off()




# FIGURE 3 -- animated gif
library(animation)
if(Sys.info()[['sysname']] == "Linux") X11.options(antialias = "none") #for Linux Systems to enable pixel-wise plotting in (animated) gif-files. 
if(Sys.info()[['sysname']] == "Windows") windows.options(antialias = "none") #for Windows Systems to enable pixel-wise plotting in (animated) gif-files. 



grazingscale <- colorRampPalette(c("black", "red2"), space = "rgb")
grazingcol <- c(grazingscale(5), color[2:3])
saveGIF( 
	for(i in seq(1, length(result$timeseries), 1/delta*2)) {
	  	par(mar = c(0,0,0,0))
		if(i < addgrazing/delta) { plot(result$timeseries[[i]], grid = FALSE, cols = color, ani = TRUE) 
		} else {
		x <- result$timeseries[[i]]
		x$cells <- as.numeric(x$cells)+3
		x$cells[x$cells == 1+3] <- (count(x, 2+3)+count(x, 3+3))[x$cells == 1+3]
		x$cells <- as.factor(x$cells+1)
		plot(x, grid = FALSE, cols = grazingcol[as.numeric(levels(x$cells))], ani = TRUE)
		}
	}
, movie.name = "herbivory_baseline.gif", img.name = "grid", convert = "convert", interval = 0.01/1,
    cmd.fun = system, clean = TRUE, ani.width = width, ani.height = height, outdir = getwd())
	