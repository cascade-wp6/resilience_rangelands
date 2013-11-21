


# specify lattice
width = 100
height = 100

states = c("+","0","-")
color <- c("black","grey80", "white") # define colors for the cell state levels

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



# time and resolution of simulation
timesteps = 5
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



res = 100
result <- data.frame(dens = seq(0,1, length=res))

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
		
	} # end of simulation.

	
parms_temp$Q_plus <- count(x_old, "+")/4
parms_temp$vul <- sum((1-parms_temp$Q_plus)[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0])/(width*height)
parms_temp$vul2 <- (1-mean(parms_temp$Q_plus[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0]))*parms_temp$rho
result$dens[i] <- parms_temp$rho
result$n[i] <- length(which(x_old$cells == "+"))
result$n_b[i] <- length(which(x_old$cells == "+"  & (1-parms_temp$Q_plus) > 0))
result$n_c[i] <- length(which(x_old$cells == "+" & (1-parms_temp$Q_plus) == 0))

result$vul[i] <- parms_temp$vul
result$vul2[i] <-  parms_temp$vul2

result$parent[i] <- with(parms_temp, (m+g/0.66))
result$stock[i] <- with(parms_temp, (m+g/rho))
result$assoc[i] <- mean(with(parms_temp, (m+g*(1-Q_plus)/0.165)[x_old$cells == "+"]))
result$assoc_b[i] <- mean(with(parms_temp, (m+g*(1-Q_plus)/0.165))[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0])
result$assoc_c[i] <- mean(with(parms_temp, (m+g*(1-Q_plus)/0.165))[x_old$cells == "+" & (1-parms_temp$Q_plus) == 0])  #vul =0.15890

result$stockassoc[i] <- mean(with(parms_temp, (m+g*(1-Q_plus)/vul))[x_old$cells == "+"])
result$stockassoc_b[i] <- mean(with(parms_temp, (m+g*(1-Q_plus)/vul))[x_old$cells == "+" & (1-parms_temp$Q_plus) > 0])
result$stockassoc_c[i] <- mean(with(parms_temp, (m+g*(1-Q_plus)/vul))[x_old$cells == "+" & (1-parms_temp$Q_plus) == 0])


}


#stopCluster(cl)

result <- result[order(result$dens),]

par(mfrow = c(1,2))

plot(NA,NA, xlim = c(0,1), ylim = c(0,10000), ylab = "intrinsic mortality rate [ind. / lattice / timestep]", xlab = "vegetation cover")

lines(parms_temp$m*n ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "death events per cell [ ind. / timestep]", xlab = "vegetation cover")

lines(parms_temp$m*n/(width*height) ~ dens, data = result, col = "black", lwd = 2)


plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "individual risk of death [ / timestep]", xlab = "vegetation cover")

lines(rep(parms_temp$m, times = length(n)) ~ dens, data = result, col = "grey50", lwd = 2)

lines(parent ~ dens, data = result, col = "black", lwd = 2)
lines(stock ~ dens, data = result, col = "red3", lwd = 2)
lines(assoc ~ dens, data = result, col = "blue", lwd = 2)
lines(assoc_b ~ dens, data = result, col = "blue", lwd = 2, lty = 2)
lines(assoc_c ~ dens, data = result, col = "blue", lwd = 2, lty = 2)
 

lines(stockassoc ~ dens, data = result, col = "green3", lwd = 2)
lines(stockassoc_b~ dens, data = result, col = "green3", lwd = 2, lty = 2)
lines(stockassoc_c~ dens, data = result, col = "green3", lwd = 2, lty = 2)



plot(NA,NA, xlim = c(0,1), ylim = c(0,10000), ylab = "individual risk of death [ / timestep]", xlab = "vegetation cover")

lines(rep(parms_temp$m, times = length(n))*n ~ dens, data = result, col = "grey50", lwd = 2)

lines(I(parent*n) ~ dens, data = result, col = "black", lwd = 2)
lines(I(stock*n) ~ dens, data = result, col = "red3", lwd = 2)
lines(I(assoc*n) ~ dens, data = result, col = "blue", lwd = 2)
lines(I(assoc_b*n) ~ dens, data = result, col = "blue", lwd = 2, lty = 2)
lines(I(assoc_c*n) ~ dens, data = result, col = "blue", lwd = 2, lty = 2)


lines((stockassoc)*n ~ dens, data = result, col = "green3", lwd = 2)
lines((stockassoc_b)*n~ dens, data = result, col = "green3", lwd = 2, lty = 2)
lines((stockassoc_c)*n~ dens, data = result, col = "green3", lwd = 2, lty = 2)


pdf("E://Dropbox//SharedWithSoniaKefi//CASCADE//13Reporting//wp6_fig.pdf", width = 9, height = 3, paper = "special")
par(mfrow = c(1,3), bg = "white", mar = c(2,2,1,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(parent*dens) ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", yaxt = "n")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(stock*dens) ~ dens, data = result, col = "black", lwd = 2)

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "deaths per timestep [ / timestep]", xlab = "vegetation cover", yaxt = "n")
lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

lines(I(stockassoc*dens) ~ dens, data = result, col = "black", lwd = 2)
lines(I(stockassoc_b*dens) ~ dens, data = result, col = "black", lwd = 2, lty = 2)
lines(I(stockassoc_c*dens) ~ dens, data = result, col = "black", lwd = 2, lty = 3)

dev.off()


pdf("E://SkyDrive//Uni//talks//talk_utrecht_fig1.pdf", width = 10, height =5, paper = "special")
par(mfrow = c(4,2), bg = "white")
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



	