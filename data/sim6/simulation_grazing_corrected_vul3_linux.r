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

#setwd("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\thread herbivory\\")
#setwd("/home/schneider/herbivory/")
#setwd("C:\\Users\\SCHNEIDER\\Documents\\CASCADE\\2013_grazingmodel\\corrected_grazing_vul\\sim6")
#setwd("E:\\Eigene Dokumente\\Uni\\projects\\2013 CASCADE grazing\\sim6")
#setwd("C:\\Users\\SCHNEIDER\\Documents\\projects\\CASCADE\\2013 CASCADE grazing\\sim6")

first_ID = 1
env <- 	unique(c(seq(0.2,1,.04)) )#seq(0.2,1,.025)
graz <- seq(.0, .2, 0.025)*0.66  #c(0.0375, 0.1125, 0.1875)
#graz <- seq(0,.1, 0.005)#seq(0,.9, 0.075)
mort <- 0.1 #seq(0.05,.15, 0.05)
init <- c("low", "high")
global <- c(TRUE, FALSE)
stock <- c(TRUE, FALSE)


lgraz <- length(graz)
lenv <- length(env)
lmort <- length(mort)
linit <- 2
lglobal <- 2
lstock <- 2
replication <- 5

#global <- rep(global, each = lgraz*lenv*lmort*linit)
#init <- rep(rep(init, each = lgraz*lenv*lmort), times = lglobal)
#mort <-	rep(rep(mort, each = lgraz*lenv), times = lglobal*linit)
#graz <-	 rep(rep(graz, each = lenv), times = lmort*lglobal*linit)
#env <-  rep(env, times = lgraz*lmort*lglobal*linit)

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


# defining parameter set for initialisation phase
parameters_ini = list(
	low = list(
		global = TRUE,
		stock = FALSE,
		starting = "low",
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
		global = TRUE,
		stock = FALSE,
		starting = "high",
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
width = 100
height = 100

# initial cell states
states = c("+","0","-")
prob = c(9/10,.9/10,0.1/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
timesteps = 3000
addgrazing = 1000
delta = 1/5

t_eval <- ((timesteps/delta)-100/delta):(timesteps/delta)+1
	
snapshots <- c(1, c(addgrazing, timesteps-20:1*50, timesteps)/delta+1)

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
			q_plus = NA,
			rho_zero = NA,
			rho_minus = NA
			)

write.table(out_header[-1,], "output.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)

header_grids <- c("ID", "global", "stock", "starting", "timestep", "g", "m0" , "b", paste("X", 1:(width*height), sep =""))

dir.create("results")

#write.table(as.data.frame(matrix(header_grids, nrow = 1, dimnames = list(1, header_grids)))[-1,], "grids.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)


########################################################################################
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

########################################################################################
########################################################################################
	
	
library(foreach)
library(doSNOW)

ubuWorker <- list(host = "kefi118", user = "schneider",
				rscript = "/usr/lib/R/bin/Rscript", #R/bin/
				snowlib = "/usr/lib/R/site-library/"
	)

#workerlist <- rep(list(ubuWorker), times = 2)
workerlist <- rep("localhost", times = 23)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)


########################################################################################
########################################################################################
	
		
foreach(iteration = parameters$ID) %dopar% {
#
#system.time( {
#iteration = 6788
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
	x_old <- initial    #old landscape for first timestep   
	
	parms_temp <- as.list(parameters[iteration,])
		# copying parms for the simulation 
	#write.table(c(iteration, i, parms_temp$g, parms_temp$m, parms_temp$b, i$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
	result$timeseries <- as.data.frame(matrix(c( parms_temp$ID, TRUE, FALSE, parms_temp$starting, 1, parms_temp$g, parms_temp$m, parms_temp$b, initial$cells), nrow = 1, dimnames = list("1", header_grids)))

	parms_temp <- parameters_ini[[c(2,1)[parameters[iteration,4]]]]
	
	result$patches <- list()
	
for(i in 2:(timesteps/delta+1)) {    #calculation loop

		if(i == addgrazing/delta ) parms_temp <- as.list(parameters[iteration,])

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 
		
		if(parms_temp$rho == 0) {flag <- FALSE} else {flag <- TRUE}
		
		
	if(flag) {
	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		
	# determine vulnerability of a random cell (summed unocupied neighbors over lattice size) 

		if(parms_temp$global == FALSE) {
				if(!parms_temp$stock) {
					parms_temp$vul <- 0.165
					} else {
					parms_temp$vul <- sum((1-parms_temp$Q_plus)[x_old$cells == "+"])/(width*height)
				}
				
				if(parms_temp$vul == 0) parms_temp$vul = 1/(width*height)
			} else {
				if(!parms_temp$stock) {
				parms_temp$vul <- 0.66
				} else {
				parms_temp$vul <- parms_temp$rho
				}
				if(parms_temp$vul == 0) parms_temp$vul = 1/(width*height)
			}
	}	else {
		parms_temp$Q_plus <- 0
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
				death <- with(parms_temp, (m+g*(1-Q_plus)/vul)*delta)
		} else {
				death <- with(parms_temp, (m+g/vul)*delta)				
		}
		death[death > 1] <- 1
	}
	
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
				
		if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
	if(flag) {		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"}
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
	if(flag) {		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"}
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
								as.data.frame(matrix(c( iteration, parms_temp$global, parms_temp$stock, parms_temp$starting, (i-1)*delta, parms_temp$g, parms_temp$m, parms_temp$b,  x_new$cells), nrow = 1, dimnames = list("1", header_grids)))
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
	
	for(i in 2:length(result$patches)) {
		result$cumpatch[paste(((snapshots-1)*delta)[-1][i])] <- sapply(bins, function(j) length(which(result$patches[[i]] >= j)) )
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
			q_plus = mean(result$q_[[1]][t_eval], na.rm = TRUE),
			rho_zero = mean(result$rho[[2]][t_eval], na.rm = TRUE),
			rho_minus = mean(result$rho[[3]][t_eval], na.rm = TRUE)
			)

write.table(out, "output.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)


#write.table(result$timeseries, "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)


gc() 
}



stopCluster(cl)


############### create pdf outputs ################


filenames <- list.files("results/")

output <- read.csv("output.csv")


pdf("patchiness_graz.pdf", height = 5, width = 6, paper = "special")
par(mfrow = c(1,1), oma = c(0,0,2,0))

with(output[output$global == FALSE,], plot(rho_plus/q_plus ~ g, pch = 20, col = highlight(b, range = c(0,1), colrange = c("black","green")), xlim = c(0,.1), ylim = c(0,1), main = "grazing at patch border" ))
 
with(output[output$global == TRUE,], plot(rho_plus/q_plus ~ g, pch = 20, col = highlight(b, range = c(0,1), colrange = c("black","green")), xlim = c(0,.1), ylim = c(0,1), main = "undifferentiated grazing" ))

dev.off()


pdf("patchiness_env.pdf", height = 5, width = 6, paper = "special")
par(mfrow = c(1,1), oma = c(0,0,2,0))

with(output[output$global == FALSE,], plot(rho_plus/q_plus ~ b, pch = 20, col = highlight(g, colrange = c("black","red2","orange","yellow"), range = c(0,.1)), xlim = c(0,1), ylim = c(0,1), main = "grazing at patch border" ))
 
with(output[output$global == TRUE,], plot(rho_plus/q_plus ~ b, pch = 20, col = highlight(g, colrange = c("black","red2","orange","yellow"), range = c(0,.1)), xlim = c(0,1), ylim = c(0,1), main = "undifferentiated grazing" ))

dev.off()


pdf("patchdist_graz.pdf", height = 9, width = 6, paper = "special")
par(mfrow = c(2,1), oma = c(0,0,2,0))

for(j in graz) {
    select <- with(parameters, ID[which(g == j & starting == "high" & global == "TRUE")])
	plot(NA, NA, xlim = c(1,10000), ylim = c(.001,1), log= "xy", main = "undifferentiated grazing" )
	color = highlight(parameters$b[parameters$ID %in% select], colrange = c("black","green3"), range = c(0,1))
	
	for(i in 1:length(select)) {
	load(paste("results/result_", select[i], sep = ""), envir = globalenv() )
	if(!is.na(result$cumpatch)) lines(I(mean/max(mean))~size, pch = 20, col = color[i], lwd = 2, data = result$cumpatch)
	}
	
	
	select <- with(parameters, ID[which(g == j & starting == "high" & global == "FALSE")])
	plot(NA, NA, xlim = c(1,10000), ylim = c(.001,1), log= "xy", main = "grazing at patch border")
	color = highlight(parameters$b[parameters$ID %in% select], colrange = c("black","green3"), range = c(0,1))
	
	for(i in 1:length(select)) {
	load(paste("results/result_", select[i], sep = ""), envir = globalenv() )
	if(!is.na(result$cumpatch)) lines(I(mean/max(mean))~size, pch = 20, col = color[i], lwd = 2, data = result$cumpatch)
	}
	
	mtext(paste("g =", j), side = 3,  outer = TRUE)
}
dev.off()


pdf("patchdist_env.pdf", height = 9, width = 6, paper = "special")
par(mfrow = c(2,1), oma = c(0,0,2,0))

for(j in  seq(.4,.6,0.02)) {
    select <- with(parameters, ID[which(b == j & starting == "high" & global == "TRUE")])
	plot(NA, NA, xlim = c(1,10000), ylim = c(.001,1), log= "xy", main = "undifferentiated grazing" )
	color = highlight(parameters$g[parameters$ID %in% select], colrange = c("black","red2","orange","yellow"), range = c(0,.1))
	
	for(i in 1:length(select)) {
	load(paste("results/result_", select[i], sep = ""), envir = globalenv() )
	if(!is.na(result$cumpatch)) lines(I(mean/max(mean))~size, pch = 20, col = color[i],lwd = 2,  data = result$cumpatch)
	}
	
	
	select <- with(parameters, ID[which(b == j & starting == "high" & global == "FALSE")])
	plot(NA, NA, xlim = c(1,10000), ylim = c(.001,1), log= "xy", main = "grazing at patch border")
	color = highlight(parameters$g[parameters$ID %in% select], colrange = c("black","red2","orange","yellow"), range = c(0,.1))
	
	for(i in 1:length(select)) {
	load(paste("results/result_", select[i], sep = ""), envir = globalenv() )
	if(!is.na(result$cumpatch)) lines(I(mean/max(mean))~size, pch = 20, col = color[i],lwd = 2,  data = result$cumpatch)
	}

	mtext(paste("b =", j), side = 3, outer = TRUE)

}
dev.off()





pdf("env_gradient.pdf", height = 6, width = 6, paper = "special")
par(mfrow = c(1,1), oma = c(0,0,2,0))
for(i in sort(unique(output$g)) ) {

	with(output[output$g == i & output$globalgrazing == TRUE,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), pch = 20, col = "black", main = paste("g = ", i)))
	
	with(output[output$g == i & output$globalgrazing == FALSE,], 
		points(b, rho_plus, pch = 20, col = "red3", cex =0.6 ))
	
}
dev.off()



pdf("grazing_gradient.pdf", height = 6, width = 6, paper = "special")
par(mfrow = c(1,1), oma = c(0,0,2,0))
for(i in sort(unique(output$b)) ) {

	with(output[output$b == i & output$globalgrazing == TRUE,], 
		plot(g, rho_plus , ylim = c(0,1), xlim = c(0,0.1), pch = 20, col = "black", main = paste("b = ", i)))
	
	with(output[output$b == i & output$globalgrazing == FALSE,], 
		points(g, rho_plus, pch = 20, col = "red3", cex = 0.6 ))
	
}
dev.off()


plotgrid_simple <- function(x , ...) { 
load(paste("results/result", x, sep = "_"))

grd <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = as.factor(states[as.integer(result$timeseries[23,-(1:7)])])  #second element contains a random row-wise, factorial vector to fill the grid 
	)
class(grd) <- c("list","landscape") # set class of object (required for plotting)

plot(grd, cols = color[match(levels(grd$cells), states)], ...)

}


plotgrid <- function(x , global, ...) { 
load(paste("results/result", x, sep = "_"))

grd <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = factor(states[as.integer(result$timeseries[23,-(1:7)])], states)  #second element contains a random row-wise, factorial vector to fill the grid 
	)
class(grd) <- c("list","landscape") # set class of object (required for plotting)

#levels(grd$cells) <- levels(grd$cells)[match( states, levels(grd$cells))]

grazingscale <- colorRampPalette(c("black", "red"), space = "rgb")
grazingcol <- c(grazingscale(32), color[2:3])


if(global) {
	vul <- length(which(grd$cells == "+"))/(width*height)
	feeding <- rep((parameters[x,]$m + parameters[x,]$g/vul), times = length(which(grd$cells == "+")))  

} else {
	vul <- sum(4-count(grd, "+" )[grd$cells == "+"])/(width*height)
	feeding <- (parameters[x,]$m + parameters[x,]$g/vul*(4-count(grd, "+" ))[grd$cells == "+"] )

}

		breaks <- seq(.1,.15, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		}


pdf("grids.pdf", width = 8, height = 8, paper = "special")
par(mfrow = c(2,2), mar = c(1,1,1,1))

plotgrid(468, TRUE)
box()

plotgrid(1350, FALSE)
box()

plotgrid(777, TRUE)
box()

plotgrid(1659, FALSE)
box()

dev.off()
