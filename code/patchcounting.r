
x <- map

#topo.colors(6, alpha = 1) rainbow(6)
rgb.palette <- colorRampPalette(c("red", "orange", "blue", "green"), space = "rgb")
cols = rgb.palette(6)[sample(1:6, 10000, replace = TRUE)]  # default color value
  
# plotting function for objects of class "landscape".
plot.map <- function(x,   dims = c(100,100), grid = FALSE, axis = FALSE, add = FALSE, ani = FALSE, ...) {
  lvls <- unique(x) 
 
  nlev <- length(lvls)
  
  if(ani & Sys.info()[['sysname']] == "Windows") adj = -0.5 else adj = 0 #this adjustment constant is added when producing a pixel accurate png or gif file, requires TRUE when the function is used to plot animated figures. 
  
  if(!add) plot(NA,NA, xlim = c(0.5+adj, dims[1]+0.5+adj), ylim = c( dims[2]+0.5+adj, 0+0.5+adj), xaxs = "i", yaxs = "i",xlab = "", ylab = "", xaxt = "n", yaxt = "n") 
 
  rect(rep(1:dims[1], times = dims[2])-.5, rep(1:dims[2], each = dims[1])-.5, rep(1:dims[1], times = dims[2])+.5, rep(1:dims[2], each = dims[1])+.5, col = cols[as.numeric(x)], border = NA)
  
 }

 
 plot(map, add = TRUE)
 plot(x_new)

 
x <- x_new
state = "+"

# get patch size and patchsize distribution
patches <- function(x, state, cumulative = TRUE) {
	pattern <- x$cells
	pattern <- pattern %in% state
	map <- rep(NA, times = prod(x$dim))
	old <- rep(99, times = prod(x$dim)) 
	
	class(map) <- c("numeric","map")
	class(old) <- c("numeric","map")
	
	while(!identical(old[pattern], map[pattern])) {
		old <- map
		count = as.integer(1)
		for(i in which(pattern)) {
			neighbors <- map[x_with_border][x_to_evaluate[i]+interact]
			if(all(is.na(neighbors)) ) { 
				map[i] <- count
			} else {
				map[i] <- min(neighbors, na.rm = TRUE)
			}
				count <- count +1
			}
 plot(map, add =TRUE)

		}
	
	map <- as.factor(map)
	patchvec <- as.vector(sapply(levels(map), function(i) length(which(map == i) ) )) 
	
	out <- vector()
	if(length(patchvec) > 0) out <- sort(patchvec) else out <- NA
	#out <- data.frame(size = unique(out), freq = sapply(unique(out), function(i) length(which(out >= i)) ))
	return(out)
	
	} 


# get patch size and patchsize distribution
patches_new <- function(x, state, cumulative = TRUE) {
	pattern <- x$cells
	pattern <- pattern %in% state
	map <- rep(NA, times = prod(x$dim))
	old <- rep(99, times = prod(x$dim)) 
	
	class(map) <- c("numeric","map")
	class(old) <- c("numeric","map")
	
	while(!identical(old[pattern], map[pattern])) {
		old <- map
		count = as.integer(1)
	
	goon <- unique(map)[!sapply(unique(map),  function(i) {
	any(sapply(which(map == i), function(k) any(map[x_with_border][x_to_evaluate[k]+interact] != i, na.rm = TRUE) ))
	}
	)]
	sapply(which(map == i), function(k) map[x_with_border][x_to_evaluate[k]+interact])
	
	
		for(i in which(pattern & map %in% goon )) {
			neighbors <- map[x_with_border][x_to_evaluate[i]+interact]
			if(all(is.na(neighbors)) ) { 
				map[i] <- count
			} else {
				map[i] <- min(neighbors, na.rm = TRUE)
			}
				count <- count +1
			}
 plot(map, add =TRUE)

		}
	
	map <- as.factor(map)
	patchvec <- as.vector(sapply(levels(map), function(i) length(which(map == i) ) )) 
	
	out <- vector()
	if(length(patchvec) > 0) out <- sort(patchvec) else out <- NA
	#out <- data.frame(size = unique(out), freq = sapply(unique(out), function(i) length(which(out >= i)) ))
	return(out)
	
	} 

	benchmark(patches(x_new, "+"), patches_new(x_new, "+"), replications = 6)
	
	patches(x_new, "+")
	
	all(patches(x_new, "+") == patches_new(x_new, "+"))