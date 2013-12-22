

rm(list=ls())


# specify lattice
width = 200
height = 200

states = c("+","0","-")
color <- c("black","grey80", "white") # define colors for the cell state levels

parms_temp <- list(
	m = 0.05,
	g = 0.03,
		b = .5 ,		# beta*eps 
		d = 0.1,		# degradation
		c_ = 0.01, 		# beta*g  
		del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
		r = 0.01, 	# regeneration rate
		f = 0.9 		# local fascilitation
)


# initial cell states
states = c("+","0","-")
prob = c(6/10,2/10,2/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

color <- c("black","grey80", "white") 
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

	

# time and resolution of simulation
timesteps = 10
delta = 1/10



count  <- function(x, neighbor, state = NULL) {
			if(is.null(state)) state = neighbor
			
			neighbors <- numeric(length = prod(x$dim))
			x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
			for(k in interact) {
				neighbors <- neighbors + x_logical_with_border[x_to_evaluate+k]
			}
			return(neighbors)	
}




res = 200
result <- data.frame(dens = seq(0,1, length=res))
x_comp <- list()

for(i in 1:res) {
j = i/res 
prob = c( i/res ,1- i/res ,0)

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
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		
		death <- with(parms_temp, m*delta)			
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		

		x_old <- x_new
		
		if(j == 0.01) x_comp[[1]] <- x_new	
		if(j == 0.05) x_comp[[2]] <- x_new		
		if(j == 0.10) x_comp[[3]] <- x_new 
		if(j == 0.30) x_comp[[4]] <- x_new 
		if(j == 0.50) x_comp[[5]] <- x_new 
		if(j == 0.60) x_comp[[6]] <- x_new 
		
	} # end of simulation.

	
parms_temp$Q_plus <- count(x_old, "+")/4

# number of cells in state +
result$n[i] <- length(which(x_old$cells == "+"))

# vulnerability_ 
parms_temp$vul <- sum((1-parms_temp$Q_plus)[x_new$cells == "+" & (1-parms_temp$Q_plus) > 0])/result$n[i]#/(width*height)

# save into result object
result$vul[i] <- parms_temp$vul

# glbal density, rho
result$dens[i] <- parms_temp$rho

# number of vegetated cells at a patch border
result$n_b[i] <- length(which(x_old$cells == "+"  & (1-parms_temp$Q_plus) > 0))

# number of vegetated cells inside a patch
result$n_c[i] <- length(which(x_old$cells == "+" & (1-parms_temp$Q_plus) == 0))
result$n_0[i] <- length(which(x_old$cells == "+" & (1-parms_temp$Q_plus) == 0))
result$n_1[i] <- length(which(x_old$cells == "+" & (1-parms_temp$Q_plus) == 1/4))
result$n_2[i] <- length(which(x_old$cells == "+" & (1-parms_temp$Q_plus) == 2/4))
result$n_3[i] <- length(which(x_old$cells == "+" & (1-parms_temp$Q_plus) == 3/4))
result$n_4[i] <- length(which(x_old$cells == "+" & (1-parms_temp$Q_plus) == 1))

# calculate average individual mortality for different models
result$parent[i] <- with(parms_temp, (g/0.5))

result$stock[i] <- with(parms_temp, (g/rho))

result$assoc3[i] <-  mean(with(parms_temp, (g/0.5*(1-Q_plus)/vul))[x_new$cells == "+"]  )  
result$assoc3_4n[i] <-  mean(with(parms_temp, (g/0.5*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 1]  ) 
result$assoc3_3n[i] <-  mean(with(parms_temp, (g/0.5*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 3/4]  ) 
result$assoc3_2n[i] <-  mean(with(parms_temp, (g/0.5*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 2/4]  ) 
result$assoc3_1n[i] <-  mean(with(parms_temp, (g/0.5*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 1/4]  )  
result$assoc3_0n[i] <-  mean(with(parms_temp, (g/0.5*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 0]  )

result$stockassoc3[i] <-  mean(with(parms_temp, (g/rho*(1-Q_plus)/vul ))[x_new$cells == "+"]  )  
result$stockassoc3_4n[i] <-  mean(with(parms_temp, (g/rho*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 4/4]  ) 
result$stockassoc3_3n[i] <-  mean(with(parms_temp, (g/rho*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 3/4]  ) 
result$stockassoc3_2n[i] <-  mean(with(parms_temp, (g/rho*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 2/4]  ) 
result$stockassoc3_1n[i] <-  mean(with(parms_temp, (g/rho*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 1/4]  )  
result$stockassoc3_0n[i] <-  mean(with(parms_temp, (g/rho*(1-Q_plus)/vul ))[x_new$cells == "+" & (1-parms_temp$Q_plus) == 0]  )

}

#stopCluster(cl)


pdf("C://Users//SCHNEIDER//SkyDrive//Uni//projects//2013 Grazing models (CASCADE)//figures//altern_models.pdf", width = 9, height = 9, paper = "special")
#pdf("E://Dropbox//SharedWithSoniaKefi//CASCADE//13Reporting//wp6_fig.pdf", width = 9, height = 3, paper = "special")
par(mfrow = c(2,2), mar = c(2,2,1,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,.5), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", xaxt = "n" )
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(parent*dens) ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,.5), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", yaxt = "n", xaxt = "n")
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(assoc3*dens) ~ dens, data = result, col = "black", lwd = 2)
lines(I(assoc3_4n*dens) ~ dens, data = result, col = "grey90", lwd = 2, lty = 2)
lines(I(assoc3_3n*dens) ~ dens, data = result, col = "grey80", lwd = 2, lty = 2)
lines(I(assoc3_2n*dens) ~ dens, data = result, col = "grey70", lwd = 2, lty = 2)
lines(I(assoc3_1n*dens) ~ dens, data = result, col = "grey60", lwd = 2, lty = 2)
lines(I(assoc3_0n*dens) ~ dens, data = result, col = "grey50", lwd = 2, lty = 3)


plot(NA,NA, xlim = c(0,1), ylim = c(0,.5), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover")
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(stock*dens) ~ dens, data = result, col = "black", lwd = 2)


plot(NA,NA, xlim = c(0,1), ylim = c(0,.5), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", yaxt = "n")
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(stockassoc3*dens) ~ dens, data = result, col = "black", lwd = 2)
lines(I(stockassoc3_4n*dens) ~ dens, data = result, col = "grey90", lwd = 2, lty = 2)
lines(I(stockassoc3_3n*dens) ~ dens, data = result, col = "grey80", lwd = 2, lty = 2)
lines(I(stockassoc3_2n*dens) ~ dens, data = result, col = "grey70", lwd = 2, lty = 2)
lines(I(stockassoc3_1n*dens) ~ dens, data = result, col = "grey60", lwd = 2, lty = 2)
lines(I(stockassoc3_0n*dens) ~ dens, data = result, col = "grey50", lwd = 2, lty = 2)



par(mfrow = c(2,2), mar = c(2,2,1,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,.5), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", xaxt = "n" )
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(parent) ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,.5), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", yaxt = "n", xaxt = "n")
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(assoc3) ~ dens, data = result, col = "black", lwd = 2)
lines(I(assoc3_4n) ~ dens, data = result, col = "grey90", lwd = 2, lty = 2)
lines(I(assoc3_3n) ~ dens, data = result, col = "grey80", lwd = 2, lty = 2)
lines(I(assoc3_2n) ~ dens, data = result, col = "grey70", lwd = 2, lty = 2)
lines(I(assoc3_1n) ~ dens, data = result, col = "grey60", lwd = 2, lty = 2)
lines(I(assoc3_0n) ~ dens, data = result, col = "grey50", lwd = 2, lty = 2)


plot(NA,NA, xlim = c(0,1), ylim = c(0,.5), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover")
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(stock) ~ dens, data = result, col = "black", lwd = 2)


plot(NA,NA, xlim = c(0,1), ylim = c(0,.5), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", yaxt = "n")
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(stockassoc3) ~ dens, data = result, col = "black", lwd = 2)
lines(I(stockassoc3_4n) ~ dens, data = result, col = "grey90", lwd = 2, lty = 2)
lines(I(stockassoc3_3n) ~ dens, data = result, col = "grey80", lwd = 2, lty = 2)
lines(I(stockassoc3_2n) ~ dens, data = result, col = "grey70", lwd = 2, lty = 2)
lines(I(stockassoc3_1n) ~ dens, data = result, col = "grey60", lwd = 2, lty = 2)
lines(I(stockassoc3_0n) ~ dens, data = result, col = "grey50", lwd = 2, lty = 2)


grazingscale <- colorRampPalette(c("black", "red"), space = "rgb")
grazingcol <- c(grazingscale(32), color[2:3])


for(s in 1:4) {
x_new <- x_comp[[s]]

par(mar = c(1,1,1,1))
grd <- x_new
	rho <- .5
	feeding <- rep((parms_temp$g/rho), times = length(which(grd$cells == "+")))  

		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		

grd <- x_new 
	Q <- count(grd, "+" )/4
	rho <- .5
	vul <- sum((1-Q)[x_new$cells == "+" & (1-Q) > 0])#/(width*height)

	feeding <- (parms_temp$g/rho*(1-Q)/vul * length(which(Q > 0)) )[grd$cells == "+"] 

		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
	
grd <- x_new
	rho <- length(which(Q > 0))/(width*height)
	feeding <- rep((parms_temp$g/rho), times = length(which(grd$cells == "+")))  

		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		


grd <- x_new 
	Q <- count(grd, "+" )/4
	rho <- length(which(Q > 0))/(width*height)
	vul <- sum((1-Q)[x_new$cells == "+" & (1-Q) > 0])#/(width*height)

	feeding <- (parms_temp$g/rho*(1-Q)/vul * length(which(Q > 0)) )[grd$cells == "+"] 

		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		
}

dev.off()





pdf("E://Dropbox//SharedWithSoniaKefi//CASCADE//13Reporting//wp6_fig.pdf", width = 9, height = 3, paper = "special")
par(mfrow = c(1,3), bg = "white", mar = c(2,2,1,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,.2), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(parent*dens) ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,.2), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", yaxt = "n")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(stock*dens) ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,.2), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", yaxt = "n")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(stockassoc*dens) ~ dens, data = result, col = "black", lwd = 2)
lines(I(stockassoc_b*dens) ~ dens, data = result, col = "black", lwd = 2, lty = 2)
lines(I(stockassoc_c*dens) ~ dens, data = result, col = "black", lwd = 2, lty = 3)

dev.off()


pdf("E://SkyDrive//Uni//talks//talk_utrecht_fig1.pdf", width = 10, height =5, paper = "special")
par(mfrow = c(1,2), bg = "white")
plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "individual risk of death [ / timestep]", xlab = "vegetation cover")

lines(rep(parms_temp$m, times = length(n)) ~ dens, data = result, col = "grey80", lwd = 2)

lines(parent ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(parent*dens) ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "individual risk of death [ / timestep]", xlab = "vegetation cover")

lines(rep(parms_temp$m, times = length(n)) ~ dens, data = result, col = "grey80", lwd = 2)
lines(parent ~ dens, data = result, col = "grey50", lwd = 2)

lines(stock ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)
lines(I(parent*dens) ~ dens, data = result, col = "grey50", lwd = 2)

lines(I(stock*dens) ~ dens, data = result, col = "black", lwd = 2)


plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "individual risk of death [ / timestep]", xlab = "vegetation cover")

lines(rep(parms_temp$m, times = length(n)) ~ dens, data = result, col = "grey80", lwd = 2)
lines(parent ~ dens, data = result, col = "grey50", lwd = 2)

lines(assoc ~ dens, data = result, col = "black", lwd = 2)
lines(assoc_b ~ dens, data = result, col = "black", lwd = 2, lty = 2)
lines(assoc_c ~ dens, data = result, col = "black", lwd = 2, lty = 2)
 

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)
lines(I(parent*dens) ~ dens, data = result, col = "grey50", lwd = 2)

lines(I(assoc*dens) ~ dens, data = result, col = "black", lwd = 2)
lines(I(assoc_b*dens) ~ dens, data = result, col = "black", lwd = 2, lty = 2)
lines(I(assoc_c*dens) ~ dens, data = result, col = "black", lwd = 2, lty = 2)


plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "individual risk of death [ / timestep]", xlab = "vegetation cover")

lines(rep(parms_temp$m, times = length(n)) ~ dens, data = result, col = "grey80", lwd = 2)
lines(parent ~ dens, data = result, col = "grey50", lwd = 2)

lines((stockassoc) ~ dens, data = result, col = "black", lwd = 2)
lines((stockassoc_b)~ dens, data = result, col = "black", lwd = 2, lty = 2)
lines((stockassoc_c)~ dens, data = result, col = "black", lwd = 2, lty = 2)


plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)
lines(I(parent*dens) ~ dens, data = result, col = "grey50", lwd = 2)

lines(I(stockassoc*dens) ~ dens, data = result, col = "black", lwd = 2)
lines(I(stockassoc_b*dens) ~ dens, data = result, col = "black", lwd = 2, lty = 2)
lines(I(stockassoc_c*dens) ~ dens, data = result, col = "black", lwd = 2, lty = 2)



dev.off()


# FIGURE 3 -- animated gif
library(animation)
if(Sys.info()[['sysname']] == "Linux") X11.options(antialias = "none") #for Linux Systems to enable pixel-wise plotting in (animated) gif-files. 
if(Sys.info()[['sysname']] == "Windows") windows.options(antialias = "none") #for Windows Systems to enable pixel-wise plotting in (animated) gif-files. 



grazingscale <- colorRampPalette(c("black", "red2"), space = "rgb")
grazingcol <- c(grazingscale(5), color[2:3])
saveGIF( 
	for(i in seq(1, length(result$timeseries), 1)) {
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
	


grazingscale <- colorRampPalette(c("black", "red"), space = "rgb")
grazingcol <- c(grazingscale(32), color[2:3])

	
pdf("E://Dropbox//SharedWithSoniaKefi//CASCADE//13Reporting//ex_grids.pdf")

par(mfrow =c(1,1), mar = c(0,0,0,0))

grd <- x_new
	vul <- .66
	feeding <- rep((0.05 + 0.033/vul), times = length(which(grd$cells == "+")))  

		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		

grd <- x_new 
	Q <- count(grd, "+" )/4
	vul <- 0.165
	feeding <- (0.05 + 0.033/vul*(1-Q)[grd$cells == "+"] )


		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
	


grd <- x_new 
	vul <- length(which(grd$cells == "+"))
	feeding <- rep((0.05 + 0.033/vul), times = length(which(grd$cells == "+")))  

		breaks <- seq(.01,.25, length = 31) #10^seq(  -1, log10(0.11 ),length = 11) #

		grd$cells <- as.numeric(grd$cells)+31
		grd$cells[grd$cells == 32] <- findInterval(feeding, breaks)+1
		grd$cells <- as.factor(grd$cells)
		levels(grd$cells) <- c(levels(grd$cells), as.character(which(! as.character(1:34) %in% levels(grd$cells))))

		plot(grd, grid = FALSE, cols = grazingcol[as.numeric(levels(grd$cells))], ani = FALSE)
		

grd <- x_new 
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



	