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

setwd("E:\\Eigene Dokumente\\Uni\\projects\\CAS01_grazing\\figures\\conceptual\\")

require(misc3d)

makeshrub <- function(id, height = 0.5, red = FALSE) {
filename <- paste("shrub",id,".png", sep = "")
if(red) filename <-paste("shrub",id,"red.png", sep = "")
png(filename, width = 325, height = 325, bg = "transparent", type =  "cairo-png" )

#bmp(paste("shrub",id,".bmp", sep = ""), width = 300, height = 300, bg = "transparent")

set.seed(id)
res = 25
x <- seq(-1,1,length=res)
z <- seq(0,2,length=res )
grd <- expand.grid(x = x, y = x, z = z)
height = 1/as.numeric(substr(id, 0,1))
width = as.numeric(substr(id, 1,2))/10
sym = 1.4
dome1 <-  1-array(grd$x^2 + grd$y^2 + grd$z^2, rep(res,3))
#struct <- array(grd$x^2 + width*grd$y^2 + height * grd$z^2, rep(length(x),3))

shapex1 <- 0.5*array( rnorm(1, mean = 2)*grd$x + rnorm(1)* grd$x^2 + rnorm(1)*grd$x^3 + rnorm(1, mean = -5)*grd$x^4, rep(res,3))
shapey1 <- 0.5*array( rnorm(1, mean = 2)*grd$y + rnorm(1)* grd$y^2 + rnorm(1)*grd$y^3 + rnorm(1, mean = -8)*grd$y^4, rep(res,3))

dome2 <-  1-array(grd$x^2 + grd$y^2 + height* grd$z^(2), rep(res,3))
#struct <- array(grd$x^2 + width*grd$y^2 + height * grd$z^2, rep(length(x),3))
#dome2[grd$z < 0] <- 0

shapex2 <- array( rnorm(1, mean = 1)*grd$x + rnorm(1)* grd$x^2 + rnorm(1)*grd$x^3 + rnorm(1, mean = -2)*grd$x^4, rep(res,3))
shapey2 <- array( rnorm(1, mean = 0.4)*grd$y + rnorm(1)* grd$y^2 + rnorm(1, mean = -2)*grd$y^3 + rnorm(1, mean = -4)*grd$y^4, rep(res,3))
#shapey2 <- array(1/sym*(rnorm(1)* -grd$y + rnorm(1)*-grd$y^2 -2*abs( rnorm(1)*-grd$y^3) ),rep(res,3))

#dome[which(dome2  * (1+shapex2) * (1+shapey2) > dome)] <- (dome2 * (1+shapex2) * (1+shapey2) )[which(dome2  * (1+shapex2) * (1+shapey2) > dome)]

noise <-  array(rnorm( length(x)*3, mean = 0, sd = 0.005), rep(length(x),3))

shrub <- pmax(dome1 + shapex1 + shapey1, dome2 + shapex2 + shapey2)+ noise

cnn <-computeContour3d(shrub, max(shrub), .75)

shrub_obj <- makeTriangles(cnn,color = "grey70", smooth = TRUE, material = "dull")

for(i in 1:3) {
shrub_obj[[i]][,1] <- shrub_obj[[i]][,1] / res - 0.5 
shrub_obj[[i]][,2] <- shrub_obj[[i]][,2] / res - 0.5 
shrub_obj[[i]][,3] <- shrub_obj[[i]][,3] / res 
}

if(red) shrub_obj <- makeTriangles(cnn,color = "red2", smooth = TRUE, material = "dull")

drawScene(shrub_obj, light = c(0.25, .5, .5,1), xlim = c(-0.5,0.5), ylim = c(-0.5,0.5), perspective = TRUE, zlim = c(0,1), col.bg = "transparent") #, xlim = c(-0.5,0.5), ylim = c(-0.5,0.5), zlim = c(-0.5,2), R.mat = cc, add = TRUE

 dev.off()

 drawScene(shrub_obj, light = c(0.25, .5, .5,1), perspective = TRUE, xlim = c(-0.5,0.5), ylim = c(-0.5,0.5), zlim = c(0,1), col.bg = "transparent") #, xlim = c(-0.5,0.5), ylim = c(-0.5,0.5), zlim = c(-0.5,2), R.mat = cc, add = TRUE

 }
 
 
 png("collection.png", width = 325*5, height = 325*5, bg = "transparent", type =  "cairo-png" )
par(mfrow = c(5,5))
 makeshrub(112)
 makeshrub(115)
 makeshrub(1164)
 makeshrub(113)
 makeshrub(164)
 

 makeshrub(212)
 makeshrub(215)
 makeshrub(2164)
 makeshrub(213)
 makeshrub(2164)
 
 makeshrub(254)
  makeshrub(255345) 
  makeshrub(2534)
  makeshrub(252)
   makeshrub(2565)
   
 makeshrub(315)
 makeshrub(3164)
 makeshrub(323)
 makeshrub(3264)
 makeshrub(344)
 
  makeshrub(345345) 
  makeshrub(3534) 
  makeshrub(3443) 
  makeshrub(3445) 
  makeshrub(354)
 
 
   dev.off()
    
 png("collection2.png", width = 325*5, height = 325*5, bg = "transparent", type =  "cairo-png" )
par(mfrow = c(5,5))
 makeshrub(14342)
 makeshrub(1454)
 makeshrub(17676)
 makeshrub(1985)
 makeshrub(1897)
 

 makeshrub(21325)
 makeshrub(213463)
 makeshrub(22124)
 makeshrub(23242)
 makeshrub(243464)
 
 makeshrub(25346)
  makeshrub(25345) 
  makeshrub(27342)
  makeshrub(28234)
   makeshrub(2964)
   
 makeshrub(3163)
 makeshrub(31342)
 makeshrub(32253)
 makeshrub(323532)
 makeshrub(3458)

 makeshrub(415)
 makeshrub(4164)
 makeshrub(413)
 makeshrub(415355)
 makeshrub(454)
 
   dev.off()
    
	
	
 png("collection3.png", width = 325*5, height = 325*5, bg = "transparent", type =  "cairo-png" )
par(mfrow = c(5,5))
 makeshrub(12642)
 makeshrub(12453)
 makeshrub(1452)
 makeshrub(2365)
 makeshrub(6341)
 

 makeshrub(214141)
 makeshrub(423)
 makeshrub(214124)
 makeshrub(234)
 makeshrub(64565)
 
 makeshrub(3453)
  makeshrub(6453) 
  makeshrub(3454)
  makeshrub(1235)
   makeshrub(42364)
   
 makeshrub(65655)
 makeshrub(64435)
 makeshrub(657453)
 makeshrub(567655)
 makeshrub(56858)

 makeshrub(51645)
 makeshrub(42464)
 makeshrub(456453)
 makeshrub(5345355)
 makeshrub(4133)
 
   dev.off()
    
	
	
cellwidth = 1/15
xpos <- seq(0,1,cellwidth) + cellwidth/2
ypos <- seq(0,1,cellwidth) + cellwidth/2
maxheight <- 0.3

shrub1 <- shrub_obj
for(i in 1:3) {
shrub1[[i]][,1] <- shrub_obj[[i]][,1] * cellwidth + xpos[1]
shrub1[[i]][,2] <- shrub_obj[[i]][,2] * cellwidth + ypos[1]
shrub1[[i]][,3] <- shrub_obj[[i]][,3] * cellwidth * maxheight



 }
 
drawScene(transformTriangles(shrub1, cc), light = c(0.25, .5, .5,1), perspective = TRUE, col.mesh = "white", zlim = c(0,1), col.bg = "transparent",  add = TRUE) #, xlim = c(-0.5,0.5), ylim = c(-0.5,0.5), zlim = c(-0.5,2)


drawScene(shrub1, light = c(0.25, .5, .5,1), perspective = TRUE, col.mesh = "white", zlim = c(0,1), col.bg = "transparent", add = TRUE) #, xlim = c(-0.5,0.5), ylim = c(-0.5,0.5), zlim = c(-0.5,2)

 
 setwd("E:\\Eigene Dokumente\\Uni\\projects\\CAS01_grazing\\data\\sim10\\")

output <- read.csv("output.csv")
#out_fits <- read.csv("output_fits.csv")
output <- output[ output$stable == TRUE, ] # & output$ID > 27541

	output$bestnum <- 0
		output$bestnum[output$best == "DES"] <- 1
		output$bestnum[output$best == "PL"] <- 2
		output$bestnum[output$best == "TPLup"] <- 3
		output$bestnum[output$best == "TPLdown"] <- 4
		output$bestnum[output$best == "EXP"] <- 5
		output$bestnum[output$best == "COV"] <- 6

excerpt <- function(x, extr = NULL) {
  out <- list()
  out$dim <- c(length(extr[1]:extr[3]),length(extr[2]:extr[4]))
  temp  <- matrix(1:length(x$cells), ncol = x$dim[2], byrow = TRUE)
  if(!is.null(extr[1]) ) temp <- temp[extr[1]:extr[3], extr[2]:extr[4], drop = FALSE]
  out$cells <-  x$cells[as.vector(temp), drop = FALSE]
  
  class(out) <- c("list", "landscape")
  
  return(out)
}



pdf("E:\\Eigene Dokumente\\Uni\\projects\\CAS01_grazing\\figures\\conceptual\\fig1_background.pdf", width = 12, height = 4, paper = "special")

# specify lattice
width = 15
height = 15

# initial cell states
states = c("+","0","-")
prob = c(6/10,2/10,2/10)
color <- c("black","grey80", "white") # define colors for the cell state levels




load(paste("results\\result", 6288, sep = "_"))
#plot(result$timeseries[[22]])
high <- excerpt(result$timeseries[[22]], extr = c(5,1,width+5,height))
pdhigh <- 	result$cumpatch[[22]]
pdhighbind <-  do.call("rbind", result$cumpatch)



load(paste("results\\result", 5462, sep = "_")) 

#plot(result$timeseries[[22]])
inter <- excerpt(result$timeseries[[22]], extr = c(45,1,width+45,height))
	pdinter <- 	result$cumpatch[[22]]
	pdinterbind <-  do.call("rbind", result$cumpatch)

load(paste("results\\result", 29878, sep = "_")) 
#plot(result$timeseries[[22]])
low <- excerpt(result$timeseries[[22]], extr = c(16,75,16+width,75+height))
	pdlow <- 	result$cumpatch[[22]]
	pdlowbind <-  do.call("rbind", result$cumpatch)
		
par(oma = c(0,0,0,0))
par(bg = "transparent")



cols = c("black","grey80", "white")

x = seq(0,1, length = width+1)
y = seq(0,1, length = height+1)
z = matrix(rep(0, times = (width+1)*(height+1)), nrow = width+1)-1.5
fill <- matrix( cols[as.numeric(high$cells)], ncol = high$dim[1])

par(new=FALSE)
persp(x, y, z, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA) -> cc

x = x+1.25
fill <- matrix( cols[as.numeric(inter$cells)], ncol = inter$dim[1])

par(new=TRUE)

persp(x, y, z, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA)

x = x+1.25
fill <- matrix( cols[as.numeric(low$cells)], ncol = low$dim[1])

par(new=TRUE)

persp(x, y, z, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA)



# 3D visualisation
x = seq(0,1, length = 2)
y = seq(0,1, length = 2)
z = matrix(rep(0, times = 4), nrow = 2)

z0 <- min(z) - .2
z <- rbind(z0, cbind(z0, z, z0), z0)
x <- c(min(x) - 1e-10, x, max(x) + 1e-10)
y <- c(min(y) - 1e-10, y, max(y) + 1e-10)

sides = "grey85"
plane = "grey92"
grids = "white"#
fill <- matrix(plane, nrow = nrow(z)-1, ncol = ncol(z)-1)
fill[c(2,8)] <- sides
fill[c(6)] <- "grey90"

par(new=TRUE)

persp(x, y, z, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = fill, scale = FALSE, border = NA, r = sqrt(2), phi = 20)

x <- x+1.25
par(new=TRUE)

persp(x, y, z, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = fill, scale = FALSE, border = NA, r = sqrt(2), phi = 20)

x <- x+1.25
par(new=TRUE)
fill[c(4)] <- "grey80"

persp(x, y, z, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = fill, scale = FALSE, border = NA, r = sqrt(2), phi = 20)
 

 
 
 
 ## helper page 1

 
cols = c("grey40","grey80", "white")

x = seq(0,1, length = width+1)
y = seq(0,1, length = height+1)
z = matrix(rep(0, times = (width+1)*(height+1)), nrow = width+1)
z2 = matrix(rep(0, times = (width+1)*(height+1)), nrow = width+1)+0.2
fill <- matrix( cols[as.numeric(high$cells)], ncol = high$dim[1])

par(new=FALSE)
persp(x, y, z, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA) -> cc
par(new=TRUE)

x = x+1.25
fill <- matrix( cols[as.numeric(inter$cells)], ncol = inter$dim[1])

par(new=TRUE)

persp(x, y, z, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA)

x = x+1.25
fill <- matrix( cols[as.numeric(low$cells)], ncol = low$dim[1])

par(new=TRUE)

persp(x, y, z, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA)

 
# page 2 
 
x = seq(0,1, length = width+1)
y = seq(0,1, length = height+1)
z = matrix(rep(0, times = (width+1)*(height+1)), nrow = width+1)
z2 = matrix(rep(0, times = (width+1)*(height+1)), nrow = width+1)+0.2
fill <- matrix( cols[as.numeric(high$cells)], ncol = high$dim[1])

par(new=FALSE)
persp(x, y, z2, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = NA, scale = FALSE, r = sqrt(2), phi = 20, border =  "grey30", lwd = 0.5) -> cc
par(new=TRUE)

x = x+1.25
fill <- matrix( cols[as.numeric(inter$cells)], ncol = inter$dim[1])

par(new=TRUE)

persp(x, y, z2, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = NA, scale = FALSE, r = sqrt(2), phi = 20, border =  "grey30", lwd = 0.5)

x = x+1.25
fill <- matrix( cols[as.numeric(low$cells)], ncol = low$dim[1])

par(new=TRUE)

persp(x, y, z2, zlim = c(-1.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,3.5), ylim = c(0,2), col = NA, scale = FALSE, r = sqrt(2), phi = 20, border = "grey30", lwd = 0.5)



dev.off()

modelcols_bg <- c("#EDE9DC","#FF0000", "#97BD54", "#C29D5E","#EE82EE","#586750") # old
modelcols_fg <- c("#BDB384","#FF0000", "#58821A", "#8B0000","#EE82EE","#30382C") # old



pdf("E:\\Eigene Dokumente\\Uni\\projects\\CAS01_grazing\\figures\\conceptual\\fig1_d.pdf", width = 7.5, height = 1.5, paper = "special", useDingbats = FALSE)

load(paste("data\\sim10\\results\\result", 5893, sep = "_")) #6288
	pdhigh <- 	result$cumpatch[[22]]
	pdhighbind <-  do.call("rbind", result$cumpatch)

load(paste("data\\sim10\\results\\result", 19720, sep = "_"))  #5462
	pdinter <- 	result$cumpatch[[22]]
	pdinterbind <-  do.call("rbind", result$cumpatch)
	
load(paste("data\\sim10\\results\\result", 5148, sep = "_"))# 29878
	pdlow <- 	result$cumpatch[[22]]
	pdlowbind <-  do.call("rbind", result$cumpatch)
		
par(mfrow = c(1,3), mar = c(3,2.2,1,1), oma = c(0,1,0,0))
plot(p~size, data = pdhighbind , log = "xy", pch = 20, xlim = c(1,10000), ylim = c(.0001,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(1, at = c(1, 100,  10000))
axis(2, at = c(.0001, .01,  1), labels =  c(.001, .01,  1) , las = 1)

PLlm <- lm(I(log(p)) ~  I(log(size)) , data = pdhighbind) 
b =  1/sum(pdhigh$n) 

TPLup <- nls(I(log(p)) ~ I( log(a) - alpha * log(size) + log(1+b/(a*size^(-alpha))) ), 
		data = pdhighbind,
		start = list(a =  exp(PLlm$coefficients[1]), alpha =  -PLlm$coefficients[2]) , #, b = 1/sum(result$cumpatch[[j]]$n)
        trace = FALSE,
		#algorithm = "port",
		#lower = c(0, 0), upper = c(1, NA),
		nls.control(maxiter = 50)
		)
		
lines(	exp(seq(log(1), log(10000), length = 100)), 
						exp(predict(TPLup, list(size = exp(seq(log(1), log(10000), length = 100 )))) ),
						col = modelcols_fg[3], lwd = 3
					)
	

plot(p~size, data = pdinterbind , log = "xy", pch = 20, xlim = c(1,10000), ylim = c(.0001,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
PLlm <- lm(I(log(p)) ~  I(log(size)) , data = pdinterbind) 
axis(1, at = c(1, 100,  10000))

TPLdown <- nls(I(log(p)) ~ I( log(a) - alpha * log(size) - (size * Sx) ), 
		data = pdinterbind,
		start = list(a = exp(PLlm$coefficients[1]), alpha =  exp(PLlm$coefficients[2]), Sx = 1/1000),
        #algorithm = "port",
		trace = FALSE
		)
		
lines(	exp(seq(log(1), log(10000), length = 100)), 
						exp(predict(TPLdown, list(size = exp(seq(log(1), log(10000), length = 100 )))) ),
						col = modelcols_fg[4], lwd = 3
					)
	

plot(p~size, data = pdlowbind , log = "xy", pch = 20, xlim = c(1,10000), ylim = c(.0001,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(1, at = c(1, 100,  10000))

dev.off()
 
 
 
 ####
 
pdf("E:\\Eigene Dokumente\\Uni\\projects\\CAS01_grazing\\figures\\conceptual\\fig2_landscape.pdf", width = 6, height = 2.5, paper = "special", useDingbats = FALSE)


 width = 10
 height = 10

load(paste("results\\result", 6288, sep = "_"))
#plot(result$timeseries[[22]])

high <- excerpt(result$timeseries[[22]], extr = c(25,60,35,70))
load(paste("results\\result", 5462, sep = "_")) 
#plot(result$timeseries[[22]])
inter <- excerpt(result$timeseries[[22]], extr = c(15,27,25,37))

par(oma = c(0,0,0,0))
par(bg = "transparent")



cols = c("black","grey80", "white")

x = seq(0,1, length = width+1)
y = seq(0,1, length = height+1)
z = matrix(rep(0, times = (width+1)*(height+1)), nrow = width+1)-0.5
fill <- matrix( cols[as.numeric(high$cells)], ncol = high$dim[1])

par(new=FALSE)
persp(x, y, z, zlim = c(-.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,2.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA) -> cc

par(new=TRUE)

x = x+1.25
fill <- matrix( cols[as.numeric(inter$cells)], ncol = inter$dim[1])
persp(x, y, z, zlim = c(-.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,2.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA) -> cc

 
# 3D visualisation page 2
x = seq(0,1, length = 2)
y = seq(0,1, length = 2)
z = matrix(rep(0, times = 4), nrow = 2)+0.3

z0 <- min(z) - .2
z <- rbind(z0, cbind(z0, z, z0), z0)
x <- c(min(x) - 1e-10, x, max(x) + 1e-10)
y <- c(min(y) - 1e-10, y, max(y) + 1e-10)

sides = "grey85"
plane = "grey92"
grids = "white"#
fill <- matrix(plane, nrow = nrow(z)-1, ncol = ncol(z)-1)
fill[c(2,8)] <- sides
fill[c(6)] <- "grey90"

par(new=FALSE)

persp(x, y, z, zlim = c(-.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,2.5), ylim = c(0,2), col = fill, scale = FALSE, border = NA, r = sqrt(2), phi = 20)

x <- x+1.25
par(new=TRUE)

persp(x, y, z, zlim = c(-.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,2.5), ylim = c(0,2), col = fill, scale = FALSE, border = NA, r = sqrt(2), phi = 20)

 
## page3
x = seq(0,1, length = width+1)
y = seq(0,1, length = height+1)
z = matrix(rep(0, times = (width+1)*(height+1)), nrow = width+1)
fill <- matrix( cols[as.numeric(high$cells)], ncol = high$dim[1])

par(new=FALSE)
persp(x, y, z, zlim = c(-.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,2.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA) -> cc

par(new=TRUE)

x = x+1.25
fill <- matrix( cols[as.numeric(inter$cells)], ncol = inter$dim[1])
persp(x, y, z, zlim = c(-.5,.5), expand = 0.4, box = FALSE, axes = FALSE, xlim = c(0,2.5), ylim = c(0,2), col = fill, scale = FALSE, r = sqrt(2), phi = 20, border = NA) -> cc



dev.off()

#####


pdf("E:\\Eigene Dokumente\\Uni\\projects\\CAS01_grazing\\figures\\conceptual\\fig2_mortalities.pdf", width =3, height = 1.5, paper = "special", useDingbats = FALSE)


layout(matrix(c(1,2,1,3,1,4), ncol = 2, byrow = TRUE), width = c(5,2))
par(mar = c(1,2,0,0), oma = c(3,1,1,3), bty = "n",  las = 1)

plot(NA,NA, xlim = c(1,0), ylim = c(0,.12), ylab = "", xlab = "", xaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.1, 2))
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)
axis(1, at  = c(0,1), labels = c("", ""), tck = 0)
## intrinsic mortality 
abline(a = 0, b = parms_temp$g/0.5, lty = 1, lwd = 3) # parent model

plot(parent ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3 ), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)
plot(parent ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3 ), ylab = "", xlab = "vegetation cover", xaxt = "n",yaxt = "n",  yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)

plot(parent ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3 ), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)



plot(NA,NA, xlim = c(1,0), ylim = c(0,.12), ylab = "", xlab = "", xaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.1, 2))
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)

abline(a = parms_temp$g, b = 0, lty = 1, lwd = 3) # livestock model
axis(1, at  = c(0,1), labels = c("", ""), tck = 0)

plot(stock ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3 ), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)
plot(stock ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)
plot(stock ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3 ), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)


plot(NA,NA, xlim = c(1,0), ylim = c(0,.12), ylab = "", xlab = "", xaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.1, 2))
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)


abline(a = 0, b = parms_temp$g/0.5, lty = 1, lwd = 3) # associative protection model
axis(1, at  = c(0,1), labels = c("", ""), tck = 0)

plot(assoc3_0n ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)
axis(4, at  = c(0,.2))
plot(assoc3_2n ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3 ), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)
plot(assoc3_4n ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)

plot(NA,NA, xlim = c(1,0), ylim = c(0,.12), ylab = "", xlab = "", xaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.1, 2))
#lines(rep(parms_temp$m, times = length(n))*dens ~ dens, data = result, col = "grey80", lwd = 2)
axis(1, at = c(0,0.5,1))
abline(a = parms_temp$g, b = 0, lty = 1, lwd = 3) # combined model

axis(1, at  = c(0,1), labels = c("", ""), tck = 0)

plot(stockassoc3_0n ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)
axis(4, at  = c(0,.2))

plot(stockassoc3_2n ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3 ), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)

plot(stockassoc3_4n ~ dens, data = result, type = "l", lwd = 3, xlim = c(1,0), ylim = c(-0.03,.3 ), ylab = "", xlab = "vegetation cover", xaxt = "n", yaxt = "n", yaxs = "i" , xaxs = "i", yaxp = c(0, 0.2, 1))
axis(4, at  = c(0,.2))
abline(h = 0)
axis(1, at  = c(0,1), labels = c("", ""), tck = 0.2)
axis(1, at  = c(0,1))

dev.off()
