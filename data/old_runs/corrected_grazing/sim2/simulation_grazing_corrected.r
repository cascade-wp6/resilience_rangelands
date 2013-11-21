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

setwd("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\thread herbivory\\corrected_grazing")

env <- 	seq(0.2,1,.025)
graz <- seq(0,.5, 0.025)
mort <- 0.1 # seq(0.1,.3, 0.1)
init <- c("low", "high")

lgraz <- length(graz)
lenv <- length(env)
lmort <- length(mort)

env <-	rep(rep(env, times = lgraz), times = lmort)
graz <-	 rep(rep(graz, each = lenv), times = lmort)
mort <-  rep(mort, each = lgraz*lenv)

replicates = 2
init <- rep(init, each = length(mort) )
env <- rep(env, times = replicates)
graz <- rep(graz, times = replicates)
mort <- rep(mort, times = replicates)

# defining parameter set
parameters = data.frame(
	ID = 1:length(graz), 
	starting = init,
	g = graz,		#grazing
	m = mort, 		# intrinsic mortality
	b = env, 		# beta*eps 
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.01, 	# regeneration rate
	f = 0.9		# local fascilitation
)


# defining parameter set
parameters_ini = list(
	low = list(
		m = 0.1, 		# intrinsic mortality
		b = 0.25, 		# beta*eps 
		d = 0.1,		# degradation
		c_ = 0.2, 		# beta*g  
		del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
		r = 0.01, 	# regeneration rate
		f = 0.9, 		# local fascilitation
		g = 0		#grazing
		), 
	high = list(
		m = 0.1, 		# intrinsic mortality
		b = 0.9, 		# beta*eps 
		d = 0.1,		# degradation
		c_ = 0.2, 		# beta*g  
		del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
		r = 0.01, 	# regeneration rate
		f = 0.9, 		# local fascilitation
		g = 0		#grazing
	)
	)


# specify lattice
width = 50
height = 50

# initial cell states
states = c("+","0","-")
prob = c(9/10,.9/10,0.1/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
timesteps = 2000
addgrazing = 500
delta = 1/5

t_eval <- ((timesteps/delta)-100/delta):(timesteps/delta)+1
	
snapshots <- c(1, c(addgrazing, timesteps-4:1*50, timesteps)/delta+1)


out_header <- data.frame(
			ID = NA,
			starting = NA, 
			grazing = NA,
			g = NA, 
			b = NA, 
			m0 = NA,
			g_gr = NA,
			m0_gr = NA,
			mortality = NA,
			mortality_border = NA,
			rho_plus = NA, 
			q_plus = NA,
			rho_zero = NA,
			rho_minus = NA
			)

write.table(out_header[-1,], "output.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)


header_grids <- c("ID", "grazing", "starting", "timestep", "g", "m" , "b", "g_gr", "m_gr", paste("X", 1:(width*height), sep =""))

dir.create("results")
dir.create("patches")



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

cols[as.integer((x-min_val)/(max_val-min_val)*steps+1)]

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

########################################################################################
########################################################################################
	
	
library(foreach)
library(doSNOW)

#	, rep(list(ubuWorker), times = 23)
workerlist <- rep("localhost", times = 11)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)

########################################################################################
########################################################################################
	

foreach(iteration = parameters$ID) %dopar% {

set.seed(iteration)
#### first run with grazing
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
			
			result$mortality <- numeric(length = timesteps/delta+1)
			#save all level global densities: rho_*
			result$rho  <- list()			
			for(i in 1:length(states)) {
				result$rho[[i]] <- numeric(length = timesteps/delta+1)
				result$rho[[i]][1]  <- sum(initial$cells == states[i])/(width*height) # write initial rho 
			}
			
			result$qplus <- numeric(length = timesteps/delta+1)
#			result$q_  <- list()			
#			for(j in 1:length(states)) {
#			result$q_[[j]] <- numeric(length = timesteps/delta)
#			result$q_[[j]][1]  <- mean(rowSums( sapply(interact, function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == "+"]+k]) )) /4# write initial local densities 
#			}
			

#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	
	parms_temp <- as.list(parameters[iteration,])
		# copying parms for the simulation 
	#write.table(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, i$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
	result$timeseries <- as.data.frame(matrix(c( parms_temp$ID, TRUE, parms_temp$starting, 1, parms_temp$g, parms_temp$m, parms_temp$b, parms_temp$g, parms_temp$m, initial$cells), nrow = 1, dimnames = list("1", header_grids)))
	
	result$patches <- list()
	
	parms_temp <- parameters_ini[[c(2,1)[parameters[iteration,2]]]]


for(i in 2:(timesteps/delta+1)) {    #calculation loop

		if(i == addgrazing/delta ) parms_temp <- as.list(parameters[iteration,])

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 

	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   

	#parms_temp$vul <- sum(parms_temp$Q_unveg[x_old$cells == "+"]*4)/(width*height)
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		death <- with(parms_temp, (m+g*(1-Q_plus))*delta)
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > .95 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		

		result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)/delta
		
		result$mortality_border[i] <- mean(death[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0], na.rm = TRUE)/delta
		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}

		result$qplus[i] <- mean(parms_temp$Q_plus[x_new$cells == "+"])
		#for(j in 1:length(states)) {
		#result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		#}
		
		if(i %in% snapshots) { 
					result$timeseries <- rbind(result$timeseries, 
								as.data.frame(matrix(c( parms_temp$ID, TRUE, parms_temp$starting, i, parms_temp$g, parms_temp$m, parms_temp$b, parms_temp$g, parms_temp$m, x_new$cells), nrow = 1, dimnames = list("1", header_grids)))
					)
								
					result$patches[[length(result$patches)+1]] <- patches(x_new,"+") 
		}
	
		
		
		x_old <- x_new
		
	} # end of simulation.

	
	# get patchsize distribution
	
	bins <- sort(unique(unlist(result$patches[-1])))
	
	
	if(length(bins) > 0) {
	result$cumpatch <- data.frame(size = bins)
	
	for(i in 2:length(result$patches)) {
		result$cumpatch[paste(((snapshots-1)*delta)[-1][i])] <- sapply(bins, function(j) length(which(result$patches[[i]] >= j)) )
	}
	
	result$cumpatch$mean <- apply(result$cumpatch[,2:length(result$patches)], 1, mean)
	result$cumpatch$sd <- apply(result$cumpatch[,2:length(result$patches)], 1, sd)
	} else {
	result$cumpatch <- NA
	}

	# save result to file
	save(result, file = paste("results\\result", iteration, "grazed", sep = "_"))
	
	#write output into csv
	out <- data.frame(
			ID = iteration,
			starting = parms_temp$starting, 
			grazing = TRUE,
			g = parms_temp$g, 
			b = parms_temp$b, 
			m0 = parms_temp$m, 	
			b0_gr = parms_temp$g,
			m0_gr = parms_temp$m,			
			mortality = mean(result$mortality[t_eval], na.rm = TRUE),
			mortality_border = mean(result$mortality_border[t_eval], na.rm = TRUE),
			rho_plus = mean(result$rho[[1]][t_eval], na.rm = TRUE), 
			#rho_plus_ungrazed = mean(result$rho[[1]][(t_eval-1)-(timesteps-addgrazing)/delta], na.rm = TRUE), 
			q_plus = mean(result$qplus[t_eval], na.rm = TRUE),
			#q_plus_ungrazed = mean(result$q_[[1]][(t_eval-1)-(timesteps-addgrazing)/delta], na.rm = TRUE),
			rho_zero = mean(result$rho[[2]][t_eval], na.rm = TRUE),
			rho_minus = mean(result$rho[[3]][t_eval], na.rm = TRUE)
			)

write.table(out, "output.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		gc()  #garbage collection




################################################################################ second run without grazing
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
			
			#result$mortality <- numeric(length = timesteps/delta+1)
			#save all level global densities: rho_*
			result$rho  <- list()			
			for(i in 1:length(states)) {
				result$rho[[i]] <- numeric(length = timesteps/delta+1)
				result$rho[[i]][1]  <- sum(initial$cells == states[i])/(width*height) # write initial rho 
			}
			
			
			result$qplus <- numeric(length = timesteps/delta+1)
			
#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	
	parms_temp <- as.list(parameters[iteration,])
		# copying parms for the simulation 
	#write.table(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, i$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
	result$timeseries <- as.data.frame(matrix(c( parms_temp$ID, FALSE, parms_temp$starting, 1, parms_temp$g, parms_temp$m, parms_temp$b, out$g, out$m0, initial$cells), nrow = 1, dimnames = list("1", header_grids)))

	result$patches <- list()
	
	parms_temp <- parameters_ini[[parameters[iteration,2]]]

for(i in 2:(timesteps/delta+1)) {    #calculation loop

		if(i == addgrazing/delta  ) {
			parms_temp <- as.list(parameters[iteration,])
			if(is.na(out$mortality)) { parms_temp$m <- out$m0 } else { parms_temp$m <- out$mortality}
			parms_temp$g <- 0
		}
		
		
		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 

	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   

	#parms_temp$vul <- sum(parms_temp$Q_unveg[x_old$cells == "+"]*4)/(width*height)
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		death <- with(parms_temp, (m)*delta)
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > .95 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		

		#result$mortality[i] <- death/delta
		
		#result$mortality_border[i] <- death/delta #mean(death[x_old$cells == "+" & parms_temp$Q_unveg > 0], na.rm = TRUE)/delta
		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}

		result$qplus[i] <- mean(parms_temp$Q_plus[x_new$cells == "+"])

		#for(j in 1:length(states)) {
		#result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		#}
		
		if(i %in% snapshots)  { 
					result$timeseries <- rbind(result$timeseries, 
								as.data.frame(matrix(c( parms_temp$ID, FALSE, parms_temp$starting, i, parms_temp$g, parms_temp$m, parms_temp$b, out$g, out$m0, x_new$cells), nrow = 1, dimnames = list("1", header_grids)))
					)
								
					result$patches[[length(result$patches)+1]] <- patches(x_new,"+") 
		}
	
		x_old <- x_new
		
	} # end of simulation.

	bins <- sort(unique(unlist(result$patches[-1])))
	
	if(length(bins) > 0) {
	result$cumpatch <- data.frame(size = bins)
	
	for(i in 2:length(result$patches)) {
		result$cumpatch[paste(((snapshots-1)*delta)[-1][i])] <- sapply(bins, function(j) length(which(result$patches[[i]] >= j)) )
	}
	
	result$cumpatch$mean <- apply(result$cumpatch[,1:length(result$patches)+1], 1, mean)
	result$cumpatch$sd <- apply(result$cumpatch[,1:length(result$patches)+1], 1, sd)
	} else {
	result$cumpatch <- NA
	}

	# save result to file
	save(result, file = paste("results\\result", iteration, "ungrazed", sep = "_"))
	

	out <- data.frame(
			ID = iteration,
			starting = parms_temp$starting, 
			grazing = FALSE,
			g = parms_temp$g, 
			b = parms_temp$b, 
			m0 = parms_temp$m, 		
			b0_gr = out$g,
			m0_gr = out$m0,						
			mortality = parms_temp$m,
			mortality_border = parms_temp$m,
			rho_plus = mean(result$rho[[1]][t_eval], na.rm = TRUE), 
			q_plus = mean(result$qplus[t_eval], na.rm = TRUE),
			rho_zero = mean(result$rho[[2]][t_eval], na.rm = TRUE),
			rho_minus = mean(result$rho[[3]][t_eval], na.rm = TRUE)
			)

write.table(out, "output.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

gc() 
}


stopCluster(cl)

filenames <- list.files("results\\")
grids <- as.data.frame(matrix(header_grids, nrow = 1, dimnames = list(1, header_grids)))[-1,]

for(i in filenames) {
load(paste("results\\", i, sep = ""))

grids <- rbind(grids, result$timeseries)

rm(result)
#file.remove(paste("results\\", i, sep = ""))
}

output <- read.csv("output.csv")


pdf("env_gradient.pdf", height = 4, width = 5, paper = "special")
par(mfrow = c(1,1), oma = c(0,0,2,0))
for(i in sort(unique(output$g)) ) {

	with(output[output$g_gr == i & output$grazing == TRUE & output$m0_gr == j,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), pch = 20, col = "red3", main = paste("g = ", i)))
	
	with(output[output$g_gr == i & output$grazing == FALSE & output$m0_gr == j,], 
		points(b, rho_plus, pch = 20, col = "black" ))
	
}
dev.off()



pdf("grazing_gradient.pdf", height = 4, width = 5, paper = "special")
par(mfrow = c(1,3), oma = c(0,0,2,0))
for(i in sort(unique(output$b)) ) {

	with(output[output$b == i & output$grazing == TRUE & output$m0_gr == j,], 
		plot(g_gr, rho_plus , ylim = c(0,1), xlim = c(0,.5), pch = 20, col = "red3", main = paste("b = ", i)))
	
	with(output[output$b == i & output$grazing == FALSE & output$m0_gr == j,], 
		points(g_gr, rho_plus, pch = 20, col = "black" ))
	
}
dev.off()




patches  <-  

for(i in filenames) {
load(paste("results\\", i, sep = ""))

patches <- rbind(grids, result$cumpatch)

rm(result)
file.remove(paste("results\\", i, sep = ""))
}


pdf("patchdist.pdf", height = 4, width = 6, paper = "special")
for(j in  unique(grids$g)[1]) {
p <- unlist( lapply( with(grids,which(grazing==1 & timestep %in% snapshots[-1:2] & g == j & m == .1 & starting == 2 & b ==  0.50)), patchdist))

plot(data.frame(size = unique(p), freq = sapply(unique(p), function(i) length(which(p == i)) )),
	log = "xy", pch = 20, xlim = c(1,2500), ylim = c(1,500),main = paste("g =", j))
}
dev.off()


