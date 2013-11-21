
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

	

# specify lattice
width = 100
height = 100

# initial cell states
states = c("+","0","-")
prob = c(9/10,.9/10,0.1/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
timesteps = 1500
delta = 1/5

t_eval <- ((timesteps/delta)-100/delta):(timesteps/delta)+1
	
snapshots <-  (timesteps-10:0*50)/delta+1

################ parameter settings
first_ID = 1  #
env <- 	seq(0.2,1,.01) #seq(0.20,0.98,.02)+.01
graz <- seq(.0, .1, 0.025) #.025 #
mort <- 0.05 #c(0.05, 0.1) #seq(0.05,.15, 0.05)
init <- round(exp(seq(log(0.25), log(0.9), length = 33))[seq(1,33,2)], digits = 3)
global <- c(FALSE, TRUE)
stock <- c(FALSE, TRUE)


lgraz <- length(graz)
lenv <- length(env)
lmort <- length(mort)
linit <- length(init)
lglobal <- length(global)
lstock <- length(stock)
replication <- 1


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





#setwd("E:\\Eigene Dokumente\\Uni\\projects\\2013 CASCADE grazing\\sim9\\")
setwd("C:\\Users\\SCHNEIDER\\Documents\\projects\\CASCADE\\2013 CASCADE grazing\\sim10\\")

filenames <- list.files("results\\")
cuts <- round( seq(0,30664, length = 11) )

library(foreach)
library(doSNOW)


#workerlist <- rep(list(ubuWorker), times = 2)
workerlist <- rep("localhost", times = 10)

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



output <- read.csv("output.csv")
output_fits <- read.csv("output_fits.csv")

output <- output[ output$g != 1 & output$stable == TRUE, ] # & output$ID > 27541

 par(mfrow = c(1,4), oma = c(0,0,2,0))
for(i in sort(unique(as.character(output$g)))[-c(1,6)] ) {
	
	with(output[output$g == 0 & output$stable == TRUE ,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), pch = 20, col = "gray80", main = paste("g = ", i)))
	
	with(output[output$g == i  & output$globalgrazing == TRUE & output$stock == FALSE & output$stable == TRUE,], 
		points(b, rho_plus , pch = 20, col = "gray20"))
	
	with(output[output$g == i & output$globalgrazing == TRUE & output$stock == TRUE & output$stable == TRUE,], 
		points(b, rho_plus, pch = 20, col = "red2", cex = 1 ))
	
	with(output[output$g == i & output$globalgrazing == FALSE & output$stock == FALSE & output$stable == TRUE,], 
		points(b, rho_plus, pch = 20, col = "darkred", cex = 1 ))
	
	with(output[output$g == i & output$globalgrazing == FALSE & output$stock == TRUE & output$stable == TRUE,], 
		points(b, rho_plus, pch = 20, col = "red", cex = 1 ))
	
}


pdf("C:\\Users\\SCHNEIDER\\SkyDrive\\Uni\\projects\\2013 Grazing models (CASCADE)\\figures\\stability.pdf", height = 7, width = 8, paper = "special")

 par(mfrow = c(2,2), mar = c(2,2,0,0), oma = c(3,3,2,1))
 
	with(output[output$g == 0 ,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), xaxt = "n", pch = 20, col = "gray80")
	)
	with(output[output$g != 0 & output$globalgrazing == TRUE & output$stock == FALSE  & output$stable == TRUE,], 
		points(b, rho_plus , pch = 20, col = highlight(g)))

	with(output[output$g == 0  & output$stable == TRUE ,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), xaxt = "n", yaxt = "n", pch = 20, col = "gray80")
	)
		#axis(4)
		
	with(output[output$g != 0 & output$globalgrazing == TRUE & output$stock == TRUE  & output$stable == TRUE,], 
		points(b, rho_plus, pch = 20, col = highlight(g), cex = 1 ))

	with(output[output$g == 0  & output$stable == TRUE ,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1),  pch = 20, col = "gray80")
	)

	with(output[output$g!= 0  & output$globalgrazing == FALSE & output$stock == FALSE  & output$stable == TRUE,], 
		points(b, rho_plus, pch = 20, col = highlight(g), cex = 1 ))
	
	with(output[output$g == 0  & output$stable == TRUE ,], 
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), yaxt = "n", pch = 20, col = "gray80")
	)
	
	with(output[output$g != 0  & output$globalgrazing == FALSE & output$stock == TRUE  & output$stable == TRUE,], 
		points(b, rho_plus, pch = 20, col = highlight(g), cex = 1 ))
		
	mtext(expression(paste("vegetation cover, ", rho["+"])), side = 2, line = 1, outer = TRUE, cex = 1.1)
	mtext(expression(paste("environmental quality, ", beta)), side = 1, line = 1, outer = TRUE, cex = 1.1)

dev.off()	




pdf("C:\\Users\\SCHNEIDER\\SkyDrive\\Uni\\projects\\2013 Grazing models (CASCADE)\\figures\\patterns.pdf", height = 7, width = 8, paper = "special")


A <- matrix(1:16, ncol = 4, byrow = TRUE)

layout(rbind( cbind(A, c(0,0,0,0), A+length(A)), rep(0, times = 9), cbind(A+length(A)*2, c(0,0,0,0), A+length(A)*3) ) , width = c(1,1,1,1,0.4,1,1,1,1), height = c(1,1,1,1,0.4,1,1,1,1))
par(oma = c(4,4,1,1)+.1, mar = c(0,0,0,0))
#layout.show(n = 48)

for(i in 1:4) {
	for(k in unique(output$g)[-1]) {
		for(j in c("0.5", "0.6", "0.7", "0.8" )) {
			#try({
				iteration = output$ID[which(output$g == k & output$b == j & output$starting %in% c(0.603, 0.653, 0.708, 0.767, 0.831, 0.900) & output$stock == c(FALSE, FALSE, TRUE, TRUE)[i] & output$globalgrazing == c(TRUE, FALSE, TRUE, FALSE)[i] & output$stable == TRUE)]
							
				load(paste("results\\result", iteration[length(iteration)], sep = "_"))
				
				#extract <- as.vector(matrix(1:10000, ncol = 100, byrow = TRUE)[1:25,1:25])
				
				x <- result$timeseries[[23]]
					
				plot(x, cols = c("grey40","grey80", "white"))
				par(new=TRUE)
				plot(NA,NA, xlim = c(1,10000), ylim = c(.00001,1), log ="xy", pch = 20, col = "#00000020", xaxt = "n", yaxt = "n", bty = "n")
					
				
						for(l in 2:22) {
							try( {
							dd3 <- data.frame(size = result$cumpatch[[l]]$size, n = result$cumpatch[[l]]$p)
							points(dd3, col = "#990000", type = "l", lwd = 2)
							}, silent = TRUE)
						}
				
				box(col = "white")
				
				par(new=FALSE)
			#})
		}
	}
}

	
dev.off()	


	
	
	

	
	
	
	
	
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

 par(mfrow = c(1,4), oma = c(0,0,2,0))

	with(output[output$g == 0.05 & output$stock == FALSE & output$globalgrazing == TRUE,], {
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), pch = 20, col = "black", main = paste("g = ", 0.075))
		arrows(b, rho_plus_ini, b, rho_plus, length = 0.05)

		}
	)
	

	with(output[output$g == 0.05 & output$stock == FALSE & output$globalgrazing == FALSE,], {
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), pch = 20, col = "black", main = paste("g = ", 0.075))
		arrows(b, rho_plus_ini, b, rho_plus, length = 0.05)

		}
	)
	
	with(output[output$g == 0.05 & output$stock == TRUE & output$globalgrazing == TRUE,], {
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), pch = 20, col = "black", main = paste("g = ", 0.075))
		arrows(b, rho_plus_ini, b, rho_plus, length = 0.05)

		}
	)
	
	with(output[output$g == 0.05 & output$stock == TRUE & output$globalgrazing == FALSE,], {
		plot(b, rho_plus , ylim = c(0,1), xlim = c(0,1), pch = 20, col = "black", main = paste("g = ", 0.075))
		arrows(b, rho_plus_ini, b, rho_plus, length = 0.05)

		}
	)
	
 
 
pdf("C:\\Users\\SCHNEIDER\\SkyDrive\\Uni\\projects\\2013 Grazing models (CASCADE)\\manuscript\\first draft\\figures\\models.pdf", height = 7, width = 12, paper = "special")

 par(mfrow = c(4,5), oma = c(2,2,2,2))

 par(mar = c(0,0,0,0) ) 
 
 for(i in 1:4) {
	for(j in sort(unique(output$g)) ){
	grazing = j
	
	temp <- output[output$g == grazing & output$stock == c(FALSE, FALSE, TRUE, TRUE)[i] & output$globalgrazing == c(TRUE, FALSE, TRUE, FALSE)[i] & output$rho_plus > 0,]
	
	 lvls <- sort(unique(temp$b))
	 repl <- sapply(lvls, function(x) length(temp$best[temp$b == x]))
	
	temp <- temp[order(temp$b), ]
	
	temp$pos1 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-(x+1)]))
	temp$pos2 <- unlist(lapply(repl, function(x) seq(0,0.01, length = x+1)[-1]))

	#temp$repl <-  repl[match(temp$b, lvls)]
	plot(NA, NA , ylim = c(0,1), xlim = c(0.2,1), bty = "n", xaxt = "n", yaxt = "n")
	rect(0.2, 0, 1, 1, col = "#F0F5CE", border = FALSE )
	
	# define colors for: desert, exp, TPL_up, PL, TPL_down, covered
	#c("#F0F5CE", "#F041D6", "#B2EB83","#FF4917", "#7171F5","#479418")
	with(temp, rect(b+pos1, 0,b+pos2,1, col = c("#F0F5CE", "violet",  "#B2EB83","red", "red3","#479418")[best+1], border = NA)
	)
	
	temp <- output[output$g == grazing & output$stock == c(FALSE, FALSE, TRUE, TRUE)[i] & output$globalgrazing == c(TRUE, FALSE, TRUE, FALSE)[i],]
	
	points(temp$b, temp$rho_plus , pch = 20, col = "black", cex = 1.25)
	#text(0.3, .9, labels = paste("g = ", grazing))
		}
}

dev.off()
	
	
 par(mfrow = c(2,2), oma = c(0,0,2,1))

 for(i in 1:4) {
	grazing = 0.05
	with(out_fits[out_fits$g == grazing & out_fits$stock == c(FALSE, FALSE, TRUE, TRUE)[i] & out_fits$globalgrazing == c(TRUE, FALSE, TRUE, FALSE)[i] & out_fits$snapshot %in% c(2000, 1900, 1800, 1700, 1600, 1500),], {
		
		plot(NA, NA , ylim = c(0,1), xlim = c(0,1), main = paste("g = ", grazing))
		
		points( b, p1 , pch = 20, col = bestmodel)
		#text(b, c(0.4,0.47,0.54, 0.6), labels = best, col = "white")
		}
	)
	}

pdf("E:\\SkyDrive\\Uni\\projects\\2013 Grazing models (CASCADE)\\figures\\alpha.pdf", height = 12, width = 7, paper = "special")
	
 par(mfrow = c(4,1), oma = c(0,0,2,1))
for(l in 1:4) {
 for(i in c(0, 0.025, 0.05, 0.075)) {
	grazing = i
	#l = 4
	with(out_fits[out_fits$g == grazing & out_fits$stock == c(FALSE, FALSE, TRUE, TRUE)[l] & out_fits$globalgrazing == c(TRUE, FALSE, TRUE, FALSE)[l] & out_fits$snapshot %in% c(2000, 1900, 1800, 1700, 1600, 1500),], {
		
		plot(NA, NA , ylim = c(0,1), xlim = c(0,1), main = paste("g = ", grazing))
		
		points( b, p1 , pch = 20, col = bestmodel)
		#text(b, c(0.4,0.47,0.54, 0.6), labels = best, col = "white")
		}
	)
	}
	}
	dev.off()

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
					try(rect(b+seq(0,0.01,length = 360)[pos-1], 0,b+seq(0,0.01,length = 360)[-pos],1, col = c("#C3F4F7",   "#F041D6", "#B2EB83","#FF4917", "#7171F5","#479418")[as.integer(best)], border = NA))

			} )
		
		}
		#arrows(b, 0,b,1, col = c("#FF000040", "#1100FF40", "#FF00FF40", "#00FFF040")[best],length = 0, lwd =8)
		
		points(b, rho_plus , pch = 20, col = "black")
		#text(b, c(0.4,0.47,0.54, 0.6), labels = best, col = "white")
		}
	)
	}

	
	
	
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
	