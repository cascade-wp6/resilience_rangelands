
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
  if(cols[1] == "auto") cols = c("black","grey80", "white")  # default color value
  
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

	
excerpt <- function(x, extr = NULL) {
  out <- list()
  out$dim <- c(length(extr[1]:extr[3]),length(extr[2]:extr[4]))
  temp  <- matrix(1:length(x$cells), ncol = x$dim[2], byrow = TRUE)
  if(!is.null(extr[1]) ) temp <- temp[extr[1]:extr[3], extr[2]:extr[4], drop = FALSE]
  out$cells <-  x$cells[as.vector(temp), drop = FALSE]
  
  class(out) <- c("list", "landscape")
  
  return(out)
}


# specify lattice
width = 100
height = 100

# initial cell states
states = c("+","0","-")
prob = c(9/10,.9/10,0.1/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
timesteps = 2000
delta = 1/5

snapshots <-  (timesteps-20:0*50)/delta+1

t_eval <- ((timesteps/delta)-100/delta):(timesteps/delta)+1
	

setwd("E:\\Eigene Dokumente\\Uni\\projects\\CAS01_grazing\\data\\sim10\\")
#setwd("C:\\Users\\SCHNEIDER\\Documents\\projects\\CAS01_grazing\\data\\sim10\\")

if(FALSE) {
filenames <- list.files("results\\")
cuts <- round( seq(0,30664, length = 11) )

library(foreach)
library(doSNOW)


#workerlist <- rep(list(ubuWorker), times = 2)
workerlist <- rep("localhost", times = 11)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)

foreach(j = 1:10, .combine = "rbind") %dopar% {
files <- filenames[(cuts[j]+1):cuts[j+1]]

output <- data.frame(
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
			rho_plus_sd = NA,
			rho_plus_ini = NA,
			q_plus = NA,
			rho_zero = NA,
			rho_minus = NA,
			largestpatch = NA, 
			largestpatch_sd = NA,
			best =  NA,
			p1 = NA,
			p2 = NA,
			p3 = NA,
			p1_sd = NA,
			p2_sd = NA,
			p3_sd = NA,
			stable = NA
	)[-1,]


	
for(i in files) {

load(paste("results\\", i, sep = ""))

output <- rbind(output, result$out)
}

return(output)
} -> output
write.table(output, "output.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)

foreach(j = 1:10, .combine = "rbind") %dopar% {
files <- filenames[(cuts[j]+1):cuts[j+1]]

out_fits <- data.frame(
		ID = NA,
		snapshot = NA,	
		starting = NA, 
		globalgrazing = NA, 
		stock = NA,
		b = NA,
		g = NA,
		m0 = NA,
		largestpatch = NA,
		bestmodel = NA,
		p1 = NA,
		p2 = NA,
		p3 = NA,		
		p1_p = NA,
		p2_p = NA,
		p3_p = NA
	)[-1,]

	
for(i in files) {

load(paste("results\\", i, sep = ""))

out_fits <- rbind(out_fits, result$fit$summary)

}

return(out_fits)
} -> out_fits
write.table(out_fits, "output_fits.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)

stopCluster(cl)

write.table(output, "output.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)
write.table(out_fits, "output_fits.csv", row.names = FALSE, col.names = TRUE, sep = ",",  append  = FALSE)


}

output <- read.csv("output.csv")
#out_fits <- read.csv("output_fits.csv")
unst_eq <- read.csv("unstable_eq.csv", sep = ",")
unst_eq <- unst_eq[unst_eq$eq >= 0.01,]
output <- output[ output$stable == TRUE, ] # & output$ID > 27541

	output$bestnum <- 0
		output$bestnum[output$best == "DES"] <- 1
		output$bestnum[output$best == "PL"] <- 2
		output$bestnum[output$best == "TPLup"] <- 3
		output$bestnum[output$best == "TPLdown"] <- 4
		output$bestnum[output$best == "EXP"] <- 5
		output$bestnum[output$best == "COV"] <- 6

		
		
# define colors for: desert, PL, TPL_up,  TPL_down, exp, covered
modelcols <- c("#FBF2BF","#FF0000", "#B2EB83", "#8B0000","#EE82EE","#479418") # old
modelcols_bg <- c("#EDE9DC","#FF0000", "#97BD54", "#C29D5E","#EE82EE","#586750") # old
modelcols_fg <- c("#BDB384","#FF0000", "#58821A", "#8B0000","#EE82EE","#30382C") # old

#modelcols <- c("#CCCCCC","#AAAAAA", "#888888", "#999999","#EE82EE","#555555") # greyscale images
#modelcols <- c("#F5FFE1","#DA77FF", "#4FC7D6", "#F7A15D","violet","#566926") #strong contrast
#modelcols <- c("#E3DFD3","#FF2C00", "#74B7F0", "#F7BC5C","violet","#586750") #weak contrast


setwd("E:\\Eigene Dokumente\\Uni\\projects\\CAS01_grazing\\")

if(FALSE) {
pdf("figures\\stability.pdf", height = 7, width = 8, paper = "special")

 par(mfrow = c(2,2), mar = c(2,2,0,0), oma = c(3,3,2,1))
 
	with(output[output$g == 0 ,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), xaxt = "n", pch = 20, col = "gray80")
	)

	with(unst_eq[unst_eq$g != 0 & unst_eq$global == TRUE & unst_eq$stock == FALSE ,], 
		points(b, eq , pch = 20, col = highlight(g,  colrange = c("#808080", "#F49D9D"))))
	with(output[output$g != 0 & output$globalgrazing == TRUE & output$stock == FALSE ,], 
		points(b, rho_plus , pch = 20, col = highlight(g)))

		
		
	with(output[output$g == 0  ,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), xaxt = "n", yaxt = "n", pch = 20, col = "gray80")
	)
	with(unst_eq[unst_eq$g != 0 & unst_eq$global == TRUE & unst_eq$stock == TRUE ,], 
		points(b, eq , pch = 20, col = highlight(g,  colrange = c("#808080", "#F49D9D"))))
	with(output[output$g != 0 & output$globalgrazing == TRUE & output$stock == TRUE  & output$stable == TRUE,], 
		points(b, rho_plus, pch = 20, col = highlight(g), cex = 1 ))

		
	with(output[output$g == 0  & output$stable == TRUE ,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1),  pch = 20, col = "gray80")
	)
	with(unst_eq[unst_eq$g != 0 & unst_eq$global == FALSE & unst_eq$stock == FALSE ,], 
		points(b, eq , pch = 20, col = highlight(g,  colrange = c("#808080", "#F49D9D"))))
	with(output[output$g!= 0  & output$globalgrazing == FALSE & output$stock == FALSE  & output$stable == TRUE,], 
		points(b, rho_plus, pch = 20, col = highlight(g), cex = 1 ))
	
	with(output[output$g == 0  & output$stable == TRUE ,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), yaxt = "n", pch = 20, col = "gray80")
	)
	with(unst_eq[unst_eq$g != 0 & unst_eq$global == FALSE & unst_eq$stock == TRUE ,], 
		points(b, eq , pch = 20, col = highlight(g,  colrange = c("#808080", "#F49D9D"))))	
	with(output[output$g != 0  & output$globalgrazing == FALSE & output$stock == TRUE  & output$stable == TRUE,], 
		points(b, rho_plus, pch = 20, col = highlight(g), cex = 1 ))
		
	mtext(expression(paste("vegetation cover, ", rho["+"])), side = 2, line = 1, outer = TRUE, cex = 1.1)
	mtext(expression(paste("environmental quality, ", beta)), side = 1, line = 1, outer = TRUE, cex = 1.1)

dev.off()	




pdf("C:\\Users\\SCHNEIDER\\Documents\\projects\\CAS01_grazing\\figures\\patterns.pdf", height =6, width = 10, paper = "special")


#A <- matrix(1:16, ncol = 4, byrow = TRUE)

#layout(rbind( cbind(A, rep(0, times = dim(A)[2]), A+length(A)), rep(0, times = 9), cbind(A+length(A)*2, c(0,0,0,0), A+length(A)*3) ) , width = c(1,1,1,1,0.4,1,1,1,1), height = c(1,1,1,1,0.4,1,1,1,1))


set.seed(232)
jvar = seq(0.5, 0.82, 0.04)+0.02
steps = length(jvar)
A <- matrix(1:(4*steps), ncol = steps, byrow = TRUE)

layout(rbind( cbind(A, rep(0, times = dim(A)[1]), A+length(A)), rep(0, times = dim(A)[2]*2+1), cbind(A+length(A)*2,  rep(0, times = dim(A)[1]), A+length(A)*3) ) , width = c(rep(1, times = steps), 0.4, rep(1, times = steps)), height = c(1,1,1,1,0.4,1,1,1,1))

par(oma = c(4,4,1,1)+.1, mar = c(0,0,0,0))
#layout.show(n = 36*4)

for(i in 1:4) {
	for(k in c(0.025, 0.05, 0.075, 0.1 )) {
		for(j in jvar) {
			#try({
				iteration = output$ID[which(output$g == k & output$b == j & output$starting %in%  sort(unique(output$starting))[12:17] & output$stock == c(FALSE, TRUE, FALSE,TRUE)[i] & output$globalgrazing == c(TRUE, TRUE, FALSE, FALSE)[i] )]
				
				powerlaw <- output$best[which(output$g == k & output$b == j & output$starting %in%  sort(unique(output$starting))[12:17] & output$stock == c(FALSE, TRUE, FALSE,TRUE)[i] & output$globalgrazing == c(TRUE, TRUE, FALSE, FALSE)[i] )] == "PL"
				
				
				#load(paste("results\\result", 5879, sep = "_"))
				
				if(length(iteration) > 1) { load(paste("data\\sim10\\results\\result", sample(iteration, 1), sep = "_")) }
				if(any(powerlaw)) { load(paste("data\\sim10\\results\\result", sample(iteration, 1, prob = powerlaw), sep = "_")) }
				if(length(iteration) == 1) { load(paste("data\\sim10\\results\\result", iteration, sep = "_")) }
				if(length(iteration) == 0) stop()
				
				#extract <- as.vector(matrix(1:10000, ncol = 100, byrow = TRUE)[1:25,1:25])
				
				#snapshot <-  which(as.character(result$fit$summary$bestmodel) == as.character(result$out$best))

				x <- excerpt(result$timeseries[[22]], extr = c(15,15,50,50))
					
				plot(x, cols = c("grey40","grey80", "white"))
				par(new=TRUE)
				plot(NA,NA, xlim = c(1,10000), ylim = c(.00001,1), log ="xy", pch = 20, col = "#00000020", xaxt = "n", yaxt = "n", bty = "n")
				
				flag_full = FALSE

				if(result$out$rho_plus < 0.02 | is.na(result$cumpatch[[22]])) {flag_desert = TRUE} else { 	# excluding deserts from model fit
									dd4 <- do.call("rbind", result$cumpatch)
									b = mean(sapply(2:22, function(x) 1/sum(result$cumpatch[[x]]$n) ))
									flag_desert = FALSE
									if(mean(sapply(2:22, function(x) length(result$cumpatch[[x]]$size))) < 3) {flag_full = TRUE} else {flag_full = FALSE}
						}
				flag = !flag_desert &  !flag_full
						
				if(flag) {
					
				
					models  <- list()
					models$AIC <- vector()
					#############
					PLlm <- lm(I(log(p)) ~  I(log(size)) , data = dd4) 
					
					
					try({models$PL <- nls(I(log(p)) ~ log(a) - alpha * log(size), 
						data = dd4,
						start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]),
						#algorithm = "port",
						trace = FALSE,
						nls.control(maxiter = 50)
						)}, silent = TRUE
					)
	
					if(!is.null(models$PL)) {
						models$AIC[1] <- AIC(models$PL)
					} else {
						models$PL  <- list(NA)
						models$AIC[1] <- NA
					}

				###########
				try({models$TPLup <- nls(I(log(p)) ~ I( log(a) - alpha * log(size) + log(1+b/(a*size^(-alpha))) ), 
						data = dd4,
						start = list(a =  exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]), #,, b = mean(sapply(2:22, function(x) 1/sum(result$cumpatch[[x]]$n) ))
  				        trace = FALSE,
						#algorithm = "port",
						#lower = c(0, 0), upper = c(1, NA),
						nls.control(maxiter = 50)
						)}, silent = TRUE
					)
	
					if(!is.null(models$TPLup)) {
						models$AIC[2] <- AIC(models$TPLup) 

					} else { 
						models$TPLup  <- list(NA)
						models$AIC[2] <- NA
					}
	
				###########
						
				try( {models$TPLdown <- nls(I(log(p)) ~ I( log(a) - alpha * log(size) - (size * Sx) ), 
						data = dd4,
						start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2], Sx = 1/1000), 
						#algorithm = "port",						
						#lower = c(0, 0, 0), upper = c(NA, NA, NA),
						trace = FALSE
						)}, silent = TRUE
					)		

					if(!is.null(models$TPLdown) & !coefficients(models$TPLdown)[3] <= 0) {
						models$AIC[3] <- AIC(models$TPLdown) 

					} else {
						models$TPLdown <- list(NA)
						models$AIC[3] <- NA
					}

				###########
	
				try( {models$EXP <- nls(I(log(p)) ~ I(log(a) -(eps*size)) , 
						data = dd4,
						start = list(a = exp(PLlm$coefficients[1]) ,eps = 1),
 				        #algorithm = "port",
						trace = FALSE
						)}, silent = TRUE
					)
	
		
					if(!is.null(models$EXP)) {
						models$AIC[4] <- AIC(models$EXP) 

					} else {
						models$EXP <- list(NA)
						models$AIC[4] <- NA
					}
	
		
					min(models$AIC, na.rm = TRUE)
					models$dAIC <- 	models$AIC -min(models$AIC, na.rm = TRUE)
	
	
					models$best <- which.min(models$AIC[-4])

					} else {
	
						if(flag_desert) { models$best <- 0 } 
						if(flag_full)	{models$best <- 5}
					}

	
					
					

				if(models$best %in% c(1,2,3,4)) {
					points(dd4$size, dd4$p, col = "#E0E43066", type = "p", lwd = 2) 
					(lines(	exp(seq(log(1), log(100000), length = 100)), 
						exp(predict(models[[models$best+1]], list(size = exp(seq(log(1), log(100000), length = 100 )))) ),
						col = modelcols[models$best+1], lwd = 4
					))
				} 
				if(models$ best == 5) {
					points(dd4$size, dd4$p, col = "#E0E43066", type = "p", lwd = 2)
					lines(c(1,mean(sapply(2:22, function(x) max(result$cumpatch[[x]]$size) ))), c(1,b ), lwd = 4, col = modelcols[models$best+1] ) 
					}
				if(models$ best == 0) {
					try({dd4 <- do.call("rbind", result$cumpatch)
					points(dd4$size, dd4$p, col = "#E0E43066", type = "p", lwd = 2)})
					}
				text(25,1, pos = 1, labels = result$out$ID, cex =1.5)
	

	

				
				box(col = "white")
				
				par(new=FALSE)
			#})
		}
	}
	

}


dev.off()	


	
	
	
	}
 
pdf("figures\\cpdf_fits.pdf", height = 7, width = 10, paper = "special", useDingbats = FALSE)

layout(matrix(c(c(1,0,0,0),2:17), ncol = 4, byrow = TRUE)) 
 par( oma = c(3,3,1,1), mar = c(2,1,1,1))

 
 temp <- output[output$g == 0 & output$rho_plus > 0,]

 
	plot(NA, NA , ylim = c(0,1), xlim = c(1,0.2), bty = "n", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
	rect(0.2, 0, 1, 1, col = modelcols_bg[1], border = FALSE )

	
	if(length(temp$ID) > 0) {
	 lvls <- sort(unique(temp$b))
	 repl <- sapply(lvls, function(x) length(temp$best[temp$b == x]))
	
	temp <- temp[order(temp$b), ]
	
	temp$pos1 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-(x+1)]))
	temp$pos2 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-1]))

	
	#temp$repl <-  repl[match(temp$b, lvls)]
	
	
	with(temp, rect(b+pos1, 0,b+pos2,1, col =  modelcols_bg[bestnum], border = NA)
	)
	}
		
	temp <- unst_eq[unst_eq$g == grazing & unst_eq$stock == c(FALSE, FALSE, TRUE, TRUE)[i] & unst_eq$global == c(TRUE, FALSE, TRUE, FALSE)[i] & unst_eq$eq > 0 & unst_eq$eq < unst_eq$rho_plus,]
	
	if(length(temp$eq) > 0) {
	#points(temp$b, temp$eq , pch = 20, col = "white", cex = 1.25)
	rect(temp$b, rep(0, length(temp$b)), c(temp$b[-1],1.01), temp$eq ,  border = NA, col = modelcols_bg[1])
	}
	
	
	temp <- output[output$g == 0,]
	
	points(temp$b, temp$rho_plus , pch = 20, col = "black", cex = 1.25)
		axis(2, at = c(0,0.5,1), las = 1)
		axis(1)	
		box()

 par(mar = c(1,1,2,1))
 
for(i in 1:4) { 
	for(j in c(0.025, 0.05, 0.075, 0.1) ){

	grazing = j
	
	temp <- output[output$g == grazing & 
	output$stock == c(FALSE, TRUE, FALSE, TRUE)[i] & 
	output$globalgrazing == c(TRUE, TRUE, FALSE, FALSE)[i] & output$rho_plus > 0,]

	plot(NA, NA , ylim = c(0,1), xlim = c(1,0.2), bty = "n", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
	rect(0.2, 0, 1, 1, col = modelcols_bg[1], border = FALSE )
	
	if(length(temp$ID) > 0) {
	 lvls <- sort(unique(temp$b))
	 repl <- sapply(lvls, function(x) length(temp$best[temp$b == x]))
	
	temp <- temp[order(temp$b), ]
	
	temp$pos1 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-(x+1)]))
	temp$pos2 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-1]))

	
	#temp$repl <-  repl[match(temp$b, lvls)]
	
	
	with(temp, rect(b+pos1, 0,b+pos2,1, col =  modelcols_bg[bestnum], border = NA)
	)
	}
	
	temp <- unst_eq[as.numeric(unst_eq$g) == grazing & unst_eq$stock == c(FALSE, TRUE, FALSE, TRUE)[i] &unst_eq$global == c(TRUE, TRUE, FALSE, FALSE)[i] & unst_eq$eq > 0 & unst_eq$eq < unst_eq$rho_plus,]

	if(length(temp$eq) > 0) {
	temp <- data.frame(b = unique(temp$b), eq = sapply(unique(temp$b), function(x) mean(temp$eq[temp$b == x]) ) )
	temp <- temp[order(temp$b),]
	
#rect(temp$b, rep(0, length(temp$b)), c(temp$b[-1],1.01), temp$eq ,  border = NA, col = "#FBF2BF")
	polygon(c(0,temp$b,1), c(0,temp$eq,0), border = NA, col = modelcols_bg[1])
	points(temp$b, temp$eq , pch = 20, col = "white", cex = 1.25)
		}
	
	
	temp <- output[output$g == grazing & output$stock == c(FALSE, TRUE, FALSE, TRUE)[i] & output$globalgrazing == c(TRUE, TRUE, FALSE, FALSE)[i],]
	
	points(temp$b, temp$rho_plus , pch = 20, col = "black", cex = 1.25)
		if(j == 0.025 ) axis(2, at = c(0,0.5,1), las = 1)	
		if(i == 4 ) axis(1)
		box()
	#text(0.3, .9, labels = paste("g = ", grazing))
		}

}

if(FALSE){
mtext("environmental quality, b", side = 1, outer = TRUE, line = 2.5,cex = 0.8)
mtext("vegetation cover, rho_plus", side = 2, outer = TRUE, line = 2.5,cex = 0.8)

mtext("type of model", side = 3, outer = TRUE, line = 2.5,cex = 1)
mtext(c("parent model", "livestock model", "associative protection", "combined"), at = c(0.1, 0.35, 0.62, 0.88), side = 3, outer = TRUE, line = 0.5,cex = 0.6)

mtext("grazing intensity, g", side = 4, outer = TRUE, line = 2.5,cex = 0.8)
mtext(c(0.1, 0.075, 0.05, 0.025), at = c(0.1, 0.35, 0.62, 0.88), side = 4, outer = TRUE, line = 0.5,cex = 0.6)
}

dev.off()
	
	
	
 
pdf("figures\\models_grids.pdf", height = 1.5, width = 7, paper = "special", useDingbats = FALSE)

	layout(matrix(c(1:5), ncol = 5, byrow = TRUE))
	
	select <- c(19586, 18440, 20449, 19720, 18906 ) # parent model
	select <- c(18367, 5492, 5067, 4654, 5545 ) # associative protecion model
	

	par( oma = c(1,1,1,1), mar = c(.5,.5,.5,.5))

 for(i in select) {
	load(paste("data\\sim10\\results\\result", i, sep = "_"))
				
				x <- excerpt(result$timeseries[[22]], extr = c(1,1,50,50))
					
				plot(x, cols = c("black","grey80", "white"))
				box()
				}

dev.off()


pdf("figures\\models_legend.pdf", height = 2.5, width = 7, paper = "special", useDingbats = FALSE)

	layout(matrix(c(1:10), ncol = 5, byrow = FALSE), height = c(1,0.66))
	
	par( oma = c(3,3,0,0))
		
 for(i in select) {
	load(paste("data\\sim10\\results\\result", i, sep = "_"))
				
				x <- excerpt(result$timeseries[[22]], extr = c(1,1,50,50))
				par(mar = c(1.5,1.5,1.5,1.5))
				plot(x, cols = c("black","grey80", "white"))
				box()
				
				par(new=FALSE, mar = c(.5,1,.5,1))
				plot(NA,NA, xlim = c(1,20000), ylim = c(.00001,1), log ="xy", bty = "l", xaxt = "n", yaxt = "n")
				axis(1, at = c(1,1000))
				if(i == select[1]) 	axis(2, at = c(.001,1), las = 1)
				flag_full = FALSE

				if(result$out$rho_plus < 0.02 | is.na(result$cumpatch[[22]])) {flag_desert = TRUE} else { 	# excluding deserts from model fit
									dd4 <- do.call("rbind", result$cumpatch)
									b = mean(sapply(2:22, function(x) 1/sum(result$cumpatch[[x]]$n) ))
									flag_desert = FALSE
									if(mean(sapply(2:22, function(x) length(result$cumpatch[[x]]$size))) < 3) {flag_full = TRUE} else {flag_full = FALSE}
						}
				flag = !flag_desert &  !flag_full
						
					models  <- list()
				if(flag) {
					
				
					models$AIC <- vector()
					#############
					PLlm <- lm(I(log(p)) ~  I(log(size)) , data = dd4) 
					
					
					try({models$PL <- nls(I(log(p)) ~ log(a) - alpha * log(size), 
						data = dd4,
						start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]),
						#algorithm = "port",
						trace = FALSE,
						nls.control(maxiter = 50)
						)}, silent = TRUE
					)
	
					if(!is.null(models$PL)) {
						models$AIC[1] <- AIC(models$PL)
					} else {
						models$PL  <- list(NA)
						models$AIC[1] <- NA
					}

				###########
				try({models$TPLup <- nls(I(log(p)) ~ I( log(a) - alpha * log(size) + log(1+b/(a*size^(-alpha))) ), 
						data = dd4,
						start = list(a =  exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]), #,, b = mean(sapply(2:22, function(x) 1/sum(result$cumpatch[[x]]$n) ))
  				        trace = FALSE,
						#algorithm = "port",
						#lower = c(0, 0), upper = c(1, NA),
						nls.control(maxiter = 50)
						)}, silent = TRUE
					)
	
					if(!is.null(models$TPLup)) {
						models$AIC[2] <- AIC(models$TPLup) 

					} else { 
						models$TPLup  <- list(NA)
						models$AIC[2] <- NA
					}
	
				###########
						
				try( {models$TPLdown <- nls(I(log(p)) ~ I( log(a) - alpha * log(size) - (size * Sx) ), 
						data = dd4,
						start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2], Sx = 1/1000), 
						#algorithm = "port",						
						#lower = c(0, 0, 0), upper = c(NA, NA, NA),
						trace = FALSE
						)}, silent = TRUE
					)		

					if(!is.null(models$TPLdown) & !coefficients(models$TPLdown)[3] <= 0) {
						models$AIC[3] <- AIC(models$TPLdown) 

					} else {
						models$TPLdown <- list(NA)
						models$AIC[3] <- NA
					}

				###########
	
				try( {models$EXP <- nls(I(log(p)) ~ I(log(a) -(eps*size)) , 
						data = dd4,
						start = list(a = exp(PLlm$coefficients[1]) ,eps = 1),
 				        #algorithm = "port",
						trace = FALSE
						)}, silent = TRUE
					)
	
		
					if(!is.null(models$EXP)) {
						models$AIC[4] <- AIC(models$EXP) 

					} else {
						models$EXP <- list(NA)
						models$AIC[4] <- NA
					}
	
		
					min(models$AIC, na.rm = TRUE)
					models$dAIC <- 	models$AIC -min(models$AIC, na.rm = TRUE)
	
	
					models$best <- which.min(models$AIC[-4])

					} else {
	
						if(flag_desert) { models$best <- 0 } 
						if(flag_full)	{models$best <- 5}
					}

	
					
					

				if(models$best %in% c(1,2,3,4)) {
					points(dd4$size, dd4$p, col = "#000000", type = "p", lwd = 2, pch = 20) 
					(lines(	exp(seq(log(1), log(100000), length = 100)), 
						exp(predict(models[[models$best+1]], list(size = exp(seq(log(1), log(100000), length = 100 )))) ),
						col = modelcols_fg[models$best+1], lwd = 4
					))
				} 
				if(models$ best == 5) {
					points(dd4$size, dd4$p, col = "#000000", type = "p", lwd = 2, pch = 20)
					lines(c(1,mean(sapply(2:22, function(x) max(result$cumpatch[[x]]$size) ))), c(1,b ), lwd = 4, col = modelcols_fg[models$best+1] ) 
					}
				
				#text(25,1, pos = 1, labels = result$out$ID, cex =1.5)
	

	
				
				par(new=FALSE)
			
		}
	
	

dev.off()
	
	
	
	
	
 
	
	
pdf("figures\\alpha.pdf", height = 5.5, width = 10, paper = "special", useDingbats = FALSE)
	
 par(mfrow= c(4,4), oma = c(3,3,2,1), mar = c(1,1,2,1))
for(l in 1:4) {
 for(i in unique(output$g)[-1]) {
	grazing = i
	#l = 4
	temp <- output[output$g == grazing & output$stock == c(FALSE, TRUE, FALSE, TRUE)[l] & output$globalgrazing == c(TRUE, TRUE, FALSE, FALSE)[l] & output$rho_plus > 0.02,]

	plot(NA, NA , ylim = c(3,0), xlim = c(1,0.2), bty = "n", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
	rect(0.2, 0, 1,3, col = modelcols_bg[1], border = FALSE )
	#if(i == 1 ) axis(2, at = c(0,-1,-2,-3))
	
	if(length(temp$ID) > 0) {
	 lvls <- sort(unique(temp$b))
	 repl <- sapply(lvls, function(x) length(temp$best[temp$b == x]))
	
	temp <- temp[order(temp$b), ]
	
	temp$pos1 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-(x+1)]))
	temp$pos2 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-1]))

	
	#temp$repl <-  repl[match(temp$b, lvls)]
	
	
	with(temp, rect(b+pos1, 0,b+pos2,3, col =  c("white", "white", "white", "white", "white", modelcols_bg[6])[bestnum], border = NA)
	)
	}
	
	temp <- output[output$g == grazing & output$stock == c(FALSE, TRUE, FALSE, TRUE)[l] & output$globalgrazing == c(TRUE, TRUE, FALSE, FALSE)[l],, drop = FALSE]

		#plot(NA, NA , ylim = c(0,3), xlim = c(0,1), main = paste("g = ", grazing))
	
		points( temp$b, temp$p1 , pch = 20, col = modelcols_fg[temp$bestnum])
		#text(b, c(0.4,0.47,0.54, 0.6), labels = best, col = "white")
		
		if(i == .025) axis(2, at = c(0,1,2,3), las = 1)	
		if(l == 4 ) axis(1)
		box()
	}
	}
	dev.off()

	
	
pdf("figures\\largestpatch.pdf", height = 5.5, width = 10, paper = "special", useDingbats = FALSE)
	
 par(mfrow= c(4,4), oma = c(3,3,2,1), mar = c(1,1,2,1))
for(l in 1:4) {
 for(i in unique(output$g)[-1]) {
	grazing = i
	#l = 4
	temp <- output[output$g == grazing & output$stock == c(FALSE, TRUE, FALSE, TRUE)[l] & output$globalgrazing == c(TRUE, TRUE, FALSE, FALSE)[l] & output$rho_plus > 0.02,]

	plot(NA, NA , ylim = c(1,10000), log = "y", xlim = c(1,0.2), bty = "n", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
	rect(0.2, 1, 1,10000, col = modelcols_bg[1], border = FALSE )
	if(i == 1 ) axis(2)
	
	if(length(temp$ID) > 0) {
	 lvls <- sort(unique(temp$b))
	 repl <- sapply(lvls, function(x) length(temp$best[temp$b == x]))
	
	temp <- temp[order(temp$b), ]
	
	temp$pos1 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-(x+1)]))
	temp$pos2 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-1]))

	
	#temp$repl <-  repl[match(temp$b, lvls)]
	
	
	with(temp, rect(b+pos1, 1,b+pos2,10000, col =  c("white", "white", "white", "white", "white", modelcols_bg[6])[bestnum], border = NA)
	)
	}
	
	
	temp <- output[output$g == grazing & output$stock == c(FALSE, TRUE, FALSE, TRUE)[l] & output$globalgrazing == c(TRUE,  TRUE, FALSE, FALSE)[l] & output$rho_plus > 0,, drop = FALSE]
		temp <- temp[order(temp$b),]
		
		#plot(NA, NA , ylim = c(1,9000), log = "y", xlim = c(0,1), xaxt = "n", yaxt = "n")
		#lower <- sapply(unique(temp$b[temp$largestpatch != 0]), function(x) min( with(temp[temp$b == x,], largestpatch - largestpatch_sd), na.rm = TRUE)  )
		#upper <- sapply(unique(temp$b[temp$largestpatch != 0]), function(x) min( with(temp[temp$b == x,], largestpatch + largestpatch_sd), na.rm = TRUE)  )
		#b <- unique(temp$b[temp$largestpatch != 0])
		
		#lines(b, lower, lty = 2, col = "grey50")
		#lines(b, upper, lty = 2, col = "grey50")

		#points( temp$b, temp$largestpatch+ temp$largestpatch_sd , pch = 20, col = paste(modelcols, "10", sep = "")[temp$bestnum])
		#points( temp$b, temp$largestpatch- temp$largestpatch_sd , pch = 20, col = paste(modelcols, "10", sep = "")[temp$bestnum])
		#text(b, c(0.4,0.47,0.54, 0.6), labels = best, col = "white")
		points( temp$b, temp$largestpatch, pch = 20, col = modelcols_fg[temp$bestnum])


		if(i == 0.025 ) axis(2, at = c(1,10,100,1000, 10000), las = 1)	
		if(l == 4 ) axis(1)
		box()
	}
	}
	dev.off()

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
 par(mfrow = c(2,1), oma = c(0,0,2,0))

for(j in  seq(.2)) {
j = 0.3
    select <- with(output, ID[which(b == j & globalgrazing == "TRUE")])
	plot(NA, NA, xlim = c(1,10000), ylim = c(.001,1), log= "xy", main = "undifferentiated grazing" )
	color = highlight(output$g[output$ID %in% select], colrange = c("black","red2","orange","yellow"), range = c(0,.1))
	
	for(i in 1:length(select)) {
	load(paste("results/result_", select[i], sep = ""), envir = globalenv() )
	if(!is.na(result$cumpatch)) lines(I(mean/max(mean))~size, pch = 20, col = color[i],lwd = 2,  data = result$cumpatch)
	}
	
	
	select <- with(output, ID[which(b == j &  globalgrazing == "FALSE")])
	plot(NA, NA, xlim = c(1,10000), ylim = c(.001,1), log= "xy", main = "grazing at patch border")
	color = highlight(output$g[output$ID %in% select], colrange = c("black","red2","orange","yellow"), range = c(0,.1))
	
	for(i in 1:length(select)) {
	load(paste("results/result_", select[i], sep = ""), envir = globalenv() )
	if(!is.na(result$cumpatch)) lines(I(mean/max(mean))~size, pch = 20, col = color[i],lwd = 2,  data = result$cumpatch)
	}

	mtext(paste("b =", j), side = 3, outer = TRUE)

}

 output[output$stock == TRUE & output$globalgrazing == TRUE & output$b == 0.925 & output$g == 0.025, ] 

 which(output$rho_plus > 0)
 #[1]  1  2  3  6  7 11 12 13 16 17

 stable <- vector()

plot(p1 ~ I(g+seq(0,0.003,length = 12)),  col = c("black", "red")[globalgrazing+1], pch = c(15, 19)[stock+1], data = out_fits)
 
 
 
 # output$ID[which(output$rho_plus > 0.4  & output$rho_plus < 0.5)] 
pdf("C:\\Users\\SCHNEIDER\\SkyDrive\\Uni\\projects\\2013 Grazing models (CASCADE)\\figures\\test.pdf", height = 12, width = 9.5, paper = "special")

 for(k in c(1,  2 , 3,  6 , 7, 11, 12 ,13, 16, 17)) {
 
#load(paste("results\\result", 38, sep = "_"))


mean(result$rho[[1]][t_eval[1:((length(t_eval)-1)/2)] ]  ) - mean(result$rho[[1]][t_eval[((length(t_eval)+1)/2):length(t_eval)] ] ) < 0.01

	dd1 <- data.frame(size = NA, n = NA )[-1,]
	for(i in 2:12 ) {
	dd1 <- rbind(dd1, data.frame(size = result$logbin[,1], n = result$logbin[,i] ))
	}
	

	dd2 <- data.frame(size = result$logbin$size, n = result$logbin$mean)



layout(matrix(c(1,2,3,3,4,4), nrow = 3, byrow = TRUE), width = c(2,1))
par(oma = c(1,1,1,1), mar = c(4,4,1,1)+.1)
plot(result$time, result$rho[[1]], ylim = c(0,1), typ = "l")

points(snapshots*delta, rep(mean(result$rho[[1]][7501:10001])+.05, times = length(snapshots)), pch = 25, col = "black" )

plot(result$timeseries[[23]], col = c("black","grey80", "white"))

#boxplot(result$rho[[1]][t_eval] ~ T)#, ylim = c(0,1))

#T = as.factor(c( rep("A", times = (length(t_eval)-1)/2), rep("B", times = (length(t_eval)+1)/2) ) )
#t.test(result$rho[[1]][t_eval] ~ T)
#summary(lm(result$rho[[1]][t_eval] ~ T))
#summary(lm(result$rho[[1]][t_eval] ~ T))$coefficients[2,4] > 0.001

#result$timeseries[!is.na(result$timeseries$ID),c(1:10)]

#x_final <- list()
#x_final$dim <- c(100, 100) 
#x_final$cells <- as.factor(result$timeseries[12,-c(1:10)]), levels = c(""))
#x_final

plot(dd1, log ="xy", pch = 20, col = "#00000020")
	points(dd2,  pch = 17, col = "black")

	
plot(NA,NA, xlim = c(1,10000), ylim = c(.00001,1), log ="xy", pch = 20, col = "#00000020")


	for(j in 2:(length(snapshots)+1) ) {
	
	#j = 3
	
	if(is.na(result$cumpatch[[j]]) || dim(result$cumpatch[[j]])[1] < 4) {flag = FALSE} else{ flag = TRUE }

	dd3 <- data.frame(size = result$cumpatch[[j]]$size, n = result$cumpatch[[j]]$p)

if(flag) {
	
	models  <- list()
	models$AIC <- vector()
	#############
PLlm <- lm(I(log(n)) ~  I(log(size)) , data = dd3) 
	
try({models$PL <- nls(I(log(n)) ~ log(a) - alpha * log(size), 
		data = dd3,
		start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]),
		trace = FALSE,
		nls.control(maxiter = 100)
		)}, silent = TRUE
	)
	
	
	if(!is.null(models$PL)) {
		models$AIC[1] <- AIC(models$PL)
	} else {
		models$PL  <- list(NA)
		models$AIC[1] <- NA
	}

###########
b=1/sum(result$cumpatch[[j]]$n)
try({models$TPLup <- nls(I(log(n)) ~ I( log(a) - alpha * log(size) + log(1+b/(a*size^(-alpha))) ), 
		data = dd3,
		start = list(a =  exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]), #, b = 1/sum(result$cumpatch[[j]]$n)
        trace = FALSE,
		#algorithm = "port",
		#lower = c(0, 0), upper = c(1, NA),
		nls.control(maxiter = 50)
		)}, silent = TRUE
	)
	
	if(!is.null(models$TPLup)) {
		models$AIC[2] <- AIC(models$TPLup) 

	} else { 
		models$TPLup  <- list(NA)
		models$AIC[2] <- NA
	}
	

###########
		
try( {models$TPLdown <- nls(I(log(n)) ~ I( log(a) - alpha * log(size) - (size/Sx) ), 
		data = dd3,
		start = list(a = exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2], Sx = 100),
        trace = FALSE
		)}, silent = TRUE
	)		

	if(!is.null(models$TPLdown)) {
		models$AIC[3] <- AIC(models$TPLdown) 

	} else {
		models$TPLdown <- list(NA)
		models$AIC[3] <- NA
	}

###########
	
try( {models$EXP <- nls(I(log(n)) ~ I(log(a) -(eps*size)) , 
		data = dd3,
		start = list(a = exp(PLlm$coefficients[1]) ,eps = 1),
        trace = FALSE
		)}, silent = TRUE
	)
	
		
	if(!is.null(models$EXP)) {
		models$AIC[4] <- AIC(models$EXP) 

	} else {
		models$EXP <- list(NA)
		models$AIC[4] <- NA
	}
	
		
	min(models$AIC, na.rm = TRUE)
	models$dAIC <- 	models$AIC -min(models$AIC, na.rm = TRUE)
	
	
	models$best <- which.min(models$AIC)

	} else {
	
		if(is.na(result$cumpatch[[j]])) { 
			models$best <- 0 
		} else { 
			models$best <- 5
		}

	
}
color = c(NA, "#FF000060", "#1100FF60", "#FF00FF60", "#FFEE0060")
if(models$best %in% c(1,2,3,4)) {
try(lines(	exp(seq(log(1), log(100000), length = 100)), 
		exp(predict(models[[models$best+1]], list(size = exp(seq(log(1), log(100000), length = 100 )))) ),
		col = color[models$best]
	))
	
	}
	points(dd3, col = "#00000060", type = "l", lwd = 2)


}	


}
dev.off()




#

			#
    ##################
##########################

##########################
##########################
##########################
##########################
### Recycle bin ##########
##########################
##########################
##########################
##########################
##########################
 ########################
 
 
 
	
 par(mfrow = c(2,5), oma = c(0,0,2,0), mar = c(3,3,2,2)+.1)
for(i in sort(unique(output$g)) ) {
	
	with(output[output$g == i & output$globalgrazing == TRUE & output$stock == FALSE,], 
		plot(rho_plus, mortality , ylim = c(0,.4), xlim = c(0,1), pch = 20, col = "black", main = paste("g = ", i)))
	
	abline(h = 0.05, col = "grey80")
	
	with(output[output$g == i & output$globalgrazing == TRUE & output$stock == TRUE,], 
		points(rho_plus, mortality, pch = 20, col = "gray50", cex = 1 ))
	
	with(output[output$g == i & output$globalgrazing == FALSE & output$stock == FALSE,], 
		points(rho_plus, mortality, pch = 20, col = "red", cex = 1 ))
	
	with(output[output$g == i & output$globalgrazing == FALSE & output$stock == TRUE,], 
		points(rho_plus, mortality, pch = 20, col = "darkred", cex = 1 ))
	
}

# par(mfrow = c(1,5), oma = c(0,0,2,0))
for(i in sort(unique(output$g)) ) {
	
	with(output[output$g == i & output$globalgrazing == TRUE & output$stock == FALSE,], 
		plot(rho_plus, mortality_border , ylim = c(0,.4), xlim = c(0,1), pch = 20, col = "black", main = paste("g = ", i)))
	
	abline(h = 0.05, col = "grey80")
	
	with(output[output$g == i & output$globalgrazing == TRUE & output$stock == TRUE,], 
		points(rho_plus, mortality_border, pch = 20, col = "gray50", cex = 1 ))
	
	with(output[output$g == i & output$globalgrazing == FALSE & output$stock == FALSE,], 
		points(rho_plus, mortality_border, pch = 20, col = "red", cex = 1 ))
	
	with(output[output$g == i & output$globalgrazing == FALSE & output$stock == TRUE,], 
		points(rho_plus, mortality_border, pch = 20, col = "darkred", cex = 1 ))
	
}

 
	#################################################
	
 par(mfrow = c(2,2), oma = c(0,0,2,1))

 for(i in 1:4) {
	grazing = 0.05
	with(output[output$g == grazing & output$stock == c(FALSE, FALSE, TRUE, TRUE)[i] & output$globalgrazing == c(TRUE, FALSE, TRUE, FALSE)[i] & output$stable == TRUE,], {
		
		plot(NA, NA , ylim = c(0,1), xlim = c(0.2,1), main = paste("g = ", grazing))
		
		for(j in ID) {
			
			with(out_fits[out_fits$ID == j, ], {
					#try(arrows(b+(0:19)/2000, 0,b+(0:19)/2000,1, col = c("#FF000010", "#1100FF10", "#FF00FF10", "#FFEE0010")[bestmodel],length = 0, lwd =.05))
					pos <- as.integer(as.factor(starting))*as.integer(as.factor(snapshot))
					#posmax = length(unique(starting))*length(unique(snapshot))
					try(rect(b+seq(0,0.01,length = 360)[pos-1], 0,b+seq(0,0.01,length = 360)[-pos],1, col = c("#FBF2BF","red", "#B2EB83", "red3","violet","#479418")[as.integer(best)], border = NA))

			} )
		
		}
		#arrows(b, 0,b,1, col = c("#FF000040", "#1100FF40", "#FF00FF40", "#00FFF040")[best],length = 0, lwd =8)
		
		points(b, rho_plus , pch = 20, col = "black")
		#text(b, c(0.4,0.47,0.54, 0.6), labels = best, col = "white")
		}
	)
	}

	
	
	