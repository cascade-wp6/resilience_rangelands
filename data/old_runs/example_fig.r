
### example
width = 20
height = 20

ex1 <- list(  
	dim = c(5, 5),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), 25, replace = T, prob = c(0.3,0.3,0.4) ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(ex1$cells) <- states  #assign cell states 
class(ex1) <- c("list","landscape") # set class of object (required for plotting)

pdf("E:\\SkyDrive\\Uni\\talks\\ex_grids.pdf")
{
par(mfrow =c(1,1), mar = c(4,4,4,4))

grd <- ex1 
	vul <- .66
	feeding <- rep((0.05 + 0.033/vul), times = length(which(grd$cells == "+")))  

		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		

grd <- ex1 
	Q <- count(grd, "+" )/4
	vul <- 0.165
	feeding <- (0.05 + 0.033/vul*(1-Q)[grd$cells == "+"] )


		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
	


grd <- ex1 
	vul <- length(which(grd$cells == "+"))
	feeding <- rep((0.05 + 0.033/vul), times = length(which(grd$cells == "+")))  

		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		

grd <- ex1 
	Q <- count(grd, "+" )/4
	vul <- sum(1-Q[grd$cells == "+" & Q > 0])/(width*height)
	feeding <- (0.05 + 0.033/vul*(1-Q)[grd$cells == "+"] )


		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
	

dev.off()
}



prob = c(9.5/10,.5/10,0/10)
prob = c(3.8/10,6.2/10,0/10)


parms_temp <- list(
	m = 0.05,
	g = 0.099,
		b = .6 ,		# beta*eps 
		d = 0.1,		# degradation
		c_ = 0.01, 		# beta*g  
		del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
		r = 0.01, 	# regeneration rate
		f = 0.9 		# local fascilitation
)


count  <- function(x, neighbor, state = NULL) {
			if(is.null(state)) state = neighbor
			
			neighbors <- numeric(length = prod(x$dim))
			x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
			for(k in interact) {
				neighbors <- neighbors + x_logical_with_border[x_to_evaluate+k]
			}
			return(neighbors)	
}


# time and resolution of simulation
timesteps = 5
delta = 1/10

grazingscale <- colorRampPalette(c("black", "red"), space = "rgb")
grazingcol <- c(grazingscale(32), color[2:3])


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

	x_old <- initial
	
for(k in 2:(timesteps/delta+1)) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- length(which(x_new$cells == "+"))/(width*height)   # get old rho plus for this timestep 

	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		parms_temp$vul <- 0.66
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		
		death <- with(parms_temp, (m+g/vul)*delta)				
		death[death > 1] <- 1		
		
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		

		x_old <- x_new
		
	} # end of simulation.

grd1 <-  x_new


	x_old <- initial
	
for(k in 2:(timesteps/delta+1)) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- length(which(x_new$cells == "+"))/(width*height)   # get old rho plus for this timestep 

	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		parms_temp$vul <- parms_temp$rho		
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		
		death <- with(parms_temp, (m+g/vul)*delta)				
		death[death > 1] <- 1		
	
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		

		x_old <- x_new
		
	} # end of simulation.

grd2 <-  x_new


	x_old <- initial
	
for(k in 2:(timesteps/delta+1)) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory	

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- length(which(x_new$cells == "+"))/(width*height)   # get old rho plus for this timestep 

	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		parms_temp$vul <- sum((1-parms_temp$Q_plus)[x_old$cells == "+"])/(width*height)
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		
		death <- with(parms_temp, (m+g*(1-Q_plus)/vul)*delta)
		death[death > 1] <- 1		
	
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		

		x_old <- x_new
		
	} # end of simulation.

grd3 <-  x_new

	
pdf("E://Dropbox//SharedWithSoniaKefi//CASCADE//13Reporting//ex_grids_high.pdf")

par(mfrow =c(1,1), mar = c(0,0,0,0))

grd <- grd1
	vul <- .66
	feeding <- rep((0.05 + 0.033/vul), times = length(which(grd$cells == "+")))  

		breaks <- seq(.010,.22, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		


grd <- grd2 
	vul <- length(which(grd$cells == "+"))/(width*height)
	feeding <- rep((0.05 + 0.033/vul), times = length(which(grd$cells == "+")))  

		breaks <- seq(.010,.22, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		

grd <- grd3
	Q <- count(grd, "+" )/4
	vul <- sum(1-Q[grd$cells == "+" & Q > 0])/(width*height)
	feeding <- (0.05 + 0.033/vul*(1-Q)[grd$cells == "+"] )


		breaks <- seq(.010,.22, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
	

dev.off()
