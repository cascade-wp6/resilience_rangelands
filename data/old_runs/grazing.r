#################################################
#
#		model: random walk
#
#
#author: Flo
#date: 19.02.2013
#
#################################################
rm(list=ls())
#setwd("/home/flo/Desktop/cellularautomata")
#setwd("E:\\Eigene Dateien\\Dropbox\\cellular automata\\")
source("E:\\Dropbox\\SharedWithMarina\\code\\cellularautomata5.r")
source("C:\\Users\\SCHNEIDER\\Dropbox\\SharedWithMarina\\code\\cellularautomata5.r")


initial <- result$timeseries[[2501]]
#initial <- landscape(100,100, levels = c("+","-","0"), prob = c(2/8,3/8))
par(mar = c(0,0,0,0))
plot(initial, grid = FALSE, col = c("black", "white", "grey80"))
#plotcounts(initial, setcountparms(initial), "+")

# time and resolution of simulation
timesteps = 100
delta = 1/10

# defining transition rules
rm(rules)
setrule("0", "+", "recolonise", rule = "(del*rho+(1-del)*Q_plus)*(b-c_*rho)", ruleset = "rules")
setrule("+", "0", "mortality", rule = "m", ruleset = "rules")
setrule("+", "0", "grazing", rule = "g*(Q_minus+Q_zero)*rho_G", ruleset = "rules")

setrule("0", "-", "degrade", rule = "d", ruleset = "rules")
setrule("-", "0", "regenerate", rule = "r + f*Q_plus", ruleset = "rules")

rules


# defining parameter set
parameters = list(
	n = 10, 		# number of goats
	g = 0.05,		# increased mortality due to individual grazing per year. goats eat 5% of body weight per day, 30 kg ~ 1.5 kg / day.
	m = 0.02, 		# intrinsic mortality
	b = 0.4, 		# beta*eps 
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.15,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.0001, 	# regeneration rate
	f = 0.9 		# local fascilitation
)

countpar <- setcountparms(initial)



# initialising result list 
result <- list()   # create an output list and
			result$dim <- initial$dim  # writes dimensions from initial object
			result$timesteps <- timesteps # writes number of intended timesteps
			result$delta <- delta # time factor for simulated steps per timestep
			result$levels <- levels(initial$cells) # writes the levels present in initial grid
			result$rules <- rules # writes the set of rules
			result$parms <- parameters # writes the initial parameters
			result$time <- seq(0, timesteps, delta) # write simulated timesteps
			
			#save all level global densities: rho_*
			result$rho  <- list()			
			for(i in 1:length(result$levels)) {
				result$rho[[i]] <- numeric(length = timesteps/delta)
				result$rho[[i]][1]  <- sum(initial$cells == result$levels[i])/prod(initial$dim) # write initial rho 
			}
			
			result$col <- c("black", "white", "grey80") # define colors for the cell state levels
			result$timeseries <- list() # create a subordinate list and
			result$timeseries[[1]] <- initial # write 'x' as the first entry of the list
			if(timeseries) for(i in 1:(timesteps/delta)+1) result$timeseries[[i]] <- initial	# allocate memory for each timeseries object

	
#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	parms_temp <- as.list(as.data.frame(parameters)) # copying parms for the simulation and multiplying
	
	if(plotting) plot(initial, grid = FALSE, col = result$col) # calls the highlevel plot


for(i in 1:(timesteps/delta)+1) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 
		parms_temp$Q_plus <- count(x_old, parms = countpar, state = "+")  # count local density of occupied fields for each cell
		parms_temp$Q_minus <- count(x_old, parms = countpar, state = "-")  # count local density of degradet fields for each cell


# 2 - drawing random numbers
		rnum <- randomnumber(prod(result$dim))  #applying a C++ function, speed gain 40%, but no assignment of seed!!!
		#rnum <- runif(prod(result$dim)) # one random number between 0 and 1 for each cell
	
	
# 3 - initialising probability boundaries for each rule and each single cell

		probs_temp <- setprobabilities(rules, parms_temp, delta )

# 4 - applying the rules to fill the cells in x_new
#  now, each rule is applied by comparing the random number for each cell to the probability boundaries:		
		for(j in 1:length(result$rules)) {
			x_new$cells[which(x_old$cells == result$rules[[j]]$is & rnum >= probs_temp$lower[[j]]& rnum < probs_temp$upper[[j]])] <- result$rules[[j]]$to
		}			
	### ------------------------------------------------------------------ ###


# 5 saving rho and local densities to result list		
		for(j in 1:length(result$levels)) {
			result$rho[[j]][i]  <- sum(x_new$cells == result$levels[j])/prod(initial$dim) # write rhovalues 
			}

# ToDo average Q_plus
	
		if(plotting) plot(x_new, add = TRUE, col = result$col) # if requested, the new state is plotted on the screen
		if(timeseries) result$timeseries[[i]] <- x_new  #if requested, the whole grid is saved to timeseries

		x_old <- x_new 
		gc()  #garbage collection
	} # end of simulation.


if(!timeseries) result$timeseries[[2]] <- x_new

class(result) <- c("list", "runresult")  # for objects of class runresult cellularautomata*.r specifies a printing method. 

result



A <- data.frame(ID = 1)
A$x <- 12
A$y <- 13
A$a <- 0
A$h <- 0

rnorm
points(A$x, A$y, pch = 21, bg = "white")



# wrapper cauchi 
library(CircStats)


for(i in 1:10) 
rwrpcauchy(1, A$a)/(2*pi)