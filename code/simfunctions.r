########################################################
# The MIT License (MIT)
#
# Copyright (c) 2014 Florian D. Schneider & Sonia Kéfi
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
########################################################


######################
## mapping function ##
######################

# required to vectorise the counting and plotting. 
# returns a map of the landscape to translate it into a vector with boundaries and another one to back-translate it to a vector without boundaries into the global environment. Needs to be called only once for the dimensions of the lattice. 


mapping <- function(width, height, boundary = "periodic", i_matrix = matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)) {
    
  # derive helper vectors for counting: 
  # transformation vector for evaluation at the border of the grid
  # set evaluation matrix 
  X <- matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
  # setting the border of the evaluation matrix X
  X <- cbind(X[,width], X, X[,1] )  
  X <- rbind(X[height,], X, X[1,] ) 
  # transformation vector which adds the border to the lattice:
  x_with_border <- as.integer(t(X))
  
  assign("x_with_border", as.integer(t(X))  , envir = .GlobalEnv )
  
  # from the matrix X (with copied border cells), which cells are the actual cells (reverse transformation vector of the previous lines) 
  #x_to_evaluate <- sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )  	
  
  
  assign("x_to_evaluate", sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )	, envir = .GlobalEnv )
  
  
  # defining the neighborhood which is to be evaluated	
  # set interaction matrix
  I <- i_matrix	
  # coordinates of neighbours in Interaction matrix I: 
  neighbours_in_I <- which(is.finite(abs(I)/abs(I)), arr.in = TRUE)
  # coordinates relative to the evaluated cell (=  which(is.na(I) ) 
  relrow <- neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
  relcol <- neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
  
  # relative position of the four direct neighbours of a cell
  #interact <- (relrow * dim(X)[2] + relcol)
  
  assign("interact", relrow * dim(X)[2] + relcol, envir = .GlobalEnv )
  
}



####################
## count function ##
####################

count  <- function(x, neighbor) {
  
  neighbors <- numeric(length = prod(x$dim))
  x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
  for(k in interact) {
    neighbors <- neighbors + x_logical_with_border[x_to_evaluate+k]
  }
  return(neighbors)  
}


#################
## patch count ##
#################

# identify and count the size of individual patches

patches <- function(x, state) {
  pattern <- x$cells
  pattern <- pattern %in% state
  map <- rep(NA, times = prod(x$dim))
  old <- rep(99, times = prod(x$dim)) 
  
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
    
  }
  
  map <- as.factor(map)
  patchvec <- as.vector(sapply(levels(map), function(i) length(which(map == i) ) )) 
  
  out <- vector()
  if(length(patchvec) > 0) out <- sort(patchvec) else out <- NA
  #out <- data.frame(size = unique(out), freq = sapply(unique(out), function(i) length(which(out >= i)) ))
  return(out)
  
} 


########################
## fitting Power Laws ##
########################


fitPL <- function(psd, p_spanning, n = NULL) {
  
  # code of fitted classes
  
  n_plants <- sum(psd$size * psd$n)/n
  
  out <- list()
  out$best <- NA
  out$AIC <- vector("numeric", length = 3)
  out$dAIC <- vector("numeric", length = 3)
  
  # criteria for vegetated state & desert state
  
  ##### linear power law model for parameter estimation
  PLlm <- lm(I(log(p)) ~  1 - I(log(size)) , data = psd) 
  
  ###########
  
  try( {out$TPLdown <- nls(I(log(p)) ~ alpha * log(size) + Sx * (1 - size) , 
                           data = psd,
                           start = list(alpha =  PLlm$coefficients, Sx = 1/1000),
                           #algorithm = "port",
                           trace = FALSE
  )}, silent = TRUE
  )    
  
  if(!is.null(out$TPLdown) & !coefficients(out$TPLdown)["Sx"] <= 0) {
    out$AIC[1] <- AIC(out$TPLdown) 
  } else {
    out$TPLdown <- list(NA)
    out$AIC[1] <- NA
  }
  
  #####
  
  try({out$PL <- nls(I(log(p)) ~ alpha * log(size), 
                     data = psd,
                     start = list( alpha =  PLlm$coefficients ),
                     trace = FALSE,
                     nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  if(!is.null(out$PL)) {
    out$AIC[2] <- AIC(out$PL)
  } else {
    out$PL  <- list(NA)
    out$AIC[2] <- NA
  }
  
  ###########
  
  try({out$TPLup <- nls(I(log(p)) ~  log(b) + log(1+(size^(alpha))/b ) , 
                        data = psd,
                        start = list( alpha =  PLlm$coefficients, b = p_spanning ) , 
                        nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  
  if(!is.null(out$TPLup)) {
    out$AIC[3] <- AIC(out$TPLup) 
  } else { 
    #result$fit$summary$TPLup  <- list(NA)
    out$TPLup  <- list(NA)
    out$AIC[3] <- NA
  }
  
  ###########
  
  out$dAIC <-   out$AIC -min(out$AIC, na.rm = TRUE)
  
  out$best <- which.min(out$AIC)+1
  
  return(out)
} 


# plotting function for objects of class "landscape".
plot.landscape <- function(x, grid = FALSE, axis = FALSE, cols = "auto", add = FALSE, ani = FALSE, ...) {
  lvls <- levels(x$cells) 
  nlev <- length(lvls)
  if(cols[1] == "auto") cols = c("black", "white")  # default color value
  
  if(ani & Sys.info()[['sysname']] == "Windows") adj = -0.5 else adj = 0 #this adjustment constant is added when producing a pixel accurate png or gif file, requires TRUE when the function is used to plot animated figures. 
  
  if(!add) plot(NA,NA, xlim = c(0.5+adj, x$dim[1]+0.5+adj), ylim = c( x$dim[2]+0.5+adj, 0+0.5+adj), bty = c("n", "o")[grid+1], xaxs = "i", yaxs = "i",xlab = "", ylab = "", xaxt = "n", yaxt = "n", ... ) 
  
  if(axis && !add) axis(3) 
  if(axis && !add) axis(2)
  
  if(grid) border = "grey80" else border = cols[as.numeric(x$cells)]
  
  rect(rep(1:x$dim[1], times = x$dim[2])-.5, rep(1:x$dim[2], each = x$dim[1])-.5, rep(1:x$dim[1], times = x$dim[2])+.5, rep(1:x$dim[2], each = x$dim[1])+.5, col = cols[as.numeric(x$cells)], border = border)
  
  if(grid) box()
}



animateCA <- function(result, filename) {
  # FIGURE 3 -- animated gif
  library(animation)
  if(Sys.info()[['sysname']] == "Linux") X11.options(antialias = "none") #for Linux Systems to enable pixel-wise plotting in (animated) gif-files. 
  if(Sys.info()[['sysname']] == "Windows") windows.options(antialias = "none") #for Windows Systems to enable pixel-wise plotting in (animated) gif-files. 
  width = result$timeseries[[1]]$dim[1]
  height = result$timeseries[[1]]$dim[2]
    
  saveGIF( 
    for(i in 1:length(result$timeseries) ) {
      par(mar = c(0,0,0,0))
      plot(result$timeseries[[i]], grid = FALSE, ani = TRUE) 
    }
    , movie.name = filename, img.name = "grid", convert = "convert", interval = 0.01/1,
    cmd.fun = system, clean = TRUE, ani.width = width, ani.height = height, outdir = getwd())
}


runCA <- function(init, parms, width = 100, height = 100, delta = 0.1, t_max = 1000, t_min = 200, t_eval = 100, isstable = 0.000001, saveeach = 50) {
  
  # cell states:
  states = c("1","0")  # vegetated, empty and degraded, respectively
  
  # mapping vectors: 
  mapping(width, height)   # create mapping vectors (see documentation)
  
  # initialize landscape: 
  parms_temp <- parms   # set temporary parameter object
  
  parms_temp$rho_one <- init   # draw initial vegetation cover
  
  init_plant <- as.integer(width*height*parms_temp$rho_one) # how many initial plants
  cells <- factor(rep("0", times = width*height), levels = states) # vector of empty cells in state "0"
  cells[sample(1:(width*height), init_plant, replace = FALSE)] <- "1"  # replace init_plant cells, randomly drawn, with "1"
  
  # wrap landscape object:
  initial <- list(  
    dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
    cells = cells  #contains a random row-wise, factorial vector to fill the grid 
  )
  levels(initial$cells) <- states  #assign cell states 
  class(initial) <- c("list","landscape") # set class of object (required for plotting)
  
  # calculate timesteps for snapshots:
  n_snaps <- length(seq(0,2*t_eval,saveeach))
  t_snaps <-  ceiling(t_max/(saveeach*n_snaps))*(n_snaps*saveeach)
  snapshots <- data.frame(time = seq(0, t_snaps, saveeach), i= seq(0, t_snaps, saveeach)/delta+1, pos = c(1:length(seq(0, t_snaps, saveeach))) )
  
  # initialise result object:
  result <- list()  # generate list object
  result$time <- seq(0, t_min, delta) # add vector of realized timesteps
  result$rho_one <- vector("numeric", length = length(result$time))  # allocate vector for global vegetation cover 
  result$rho_one[1] <- parms_temp$rho_one   # fill in first value of vegetation cover
  result$q_one_one  <- vector("numeric", length = length(result$time))  # allocate vector for average local cover
  result$q_one_one[1]  <- mean(subset(count(initial, "1")/4, initial$cells == "1") )  # fill in first vector of average local cover
  
  result$timeseries <- list()
  result$timeseries[[1]] <- initial
  
  # --------------- simulation -----------------
  
  # initialise simulation variables: 
  x_old <- initial  # ghost matrix at t_i
  x_new <- x_old
  stability <- 1  # check value for stability
  i = 0  # iterator for simulation timesteps
  
  # starting iterations:
  while(stability > isstable & i <= t_max/delta & parms_temp$rho_one > 0) {
    
    i <- i +1  # increase iterator
   
    # model specific part:
    # 1 - setting time-step parameters
    parms_temp$q_one_one <- count(x_old, "1")/4  # count local density of occupied fields for each cell:
    
    # 2 - drawing random numbers
    rnum <- runif(width*height) # one random number between 0 and 1 for each cell
    
    # 3 - setting transition probabilities
    growth <- with(parms_temp, (r * (b + (1-b)*f*q_one_one) * rho_one^(1 + alpha) * ( 1 - (rho_one / (K * (1-c*q_one_one) ))) / (1 - rho_one))  *delta)  # recolonisation rates of all cells 
    
    growth[growth < 0] <- 0
    
    death <- with(parms_temp, (m + ( (a + v*q_one_one) * L * (1-p*q_one_one) * rho_one^(q) )/( 1 + (a + v*q_one_one) * h * rho_one^(1+q) )) *delta)   # set probability of death for each cell
    
    death[death < 0] <- 0
    
    # check for sum of probabilities to be inferior 1 and superior 0
    if(any(c(growth, death) > 1 )) warning(paste("a set probability is exceeding 1 in time step", i, "! decrease delta!!!")) 
    #if(any(c(growth, death) < 0)) warning(paste("a set probability falls below 0 in time step", i, "! balance parameters!!!")) 
    
    # 4 - apply transition probabilities  
    
    x_new$cells[which(x_old$cells == "0" & rnum <= growth)] <- "1"
    x_new$cells[which(x_old$cells == "1" & rnum <= death)] <- "0"
    
    # 5 - saving state of the new grid		
    
    parms_temp$rho_one <- sum(x_new$cells == "1")/(width*height) # store current global vegetation cover
    result$rho_one[i] <- parms_temp$rho_one # write current vegetation cover to result
    #result$q_one_one[i]  <- mean(subset(count(x_new, "1")/4, x_new$cells == "1") ) # write current average local cover to result
    x_old <- x_new # replace ghost matrix for next iteration
    
    
    if(i %in% snapshots$i) {  
      result$timeseries[[match(i, snapshots$i)]] <- x_new
    }
    
    
    if(i > t_min/delta+1) { # if we are over the minimal timespan 
      t_1 <- (i-2*t_eval/delta):(i-t_eval/delta)-1 # vector of t_eval timesteps previous to the last t_eval timesteps
      t_2 <- (i-t_eval/delta):(i) # vector of the last t_eval timesteps 
      
      if(parms_temp$rho_one > 0) { 
        stability <- (abs(mean(result$rho_one[t_1]) - mean(result$rho_one[t_2])))/(mean(result$rho_one[t_1])) # calculate stability, i.e. difference between the mean cover in the two evaluation periods t_1 and t_2 
      } else {
        stability <- 0 # set stability to 0 if cover is 0, immediate stop of simulation
      }
      
      result$time[i] <- i*delta # save timestep to results
      
    }
    
  } 
  
  return(result)
}



### Ordinary differential equations

## global model

d_rho <- function(rho, parms, z = 4) { 
  with(parms, r* b* rho * (1 - rho/K) - m * rho - (a * rho * L)/(1 + a * h * rho) )  
}


odesys <- function (t, rho, parms = model_parms) {
  list( d_rho(rho, parms)  )
}

## mean-field


d_rho_mean <- function(rho, parms) { 
  growth <- with(parms, r *  rho^( 1 + alpha) * (b + (1-b) * f * rho)  * (1 - rho/(K * (1-c*rho) ) ))
  if(growth <= 0| is.na(growth)) {growth <- 0}
  
  mortality <- with(parms,  m * rho + ((a+v*rho) * rho^( 1 + q) * L * (1-p*rho))/(1 + (a+v*rho) * h * rho^( 1 + q)) )  
  if(mortality <= 0 | is.na(mortality)) {mortality <- 0}
  
  return(growth - mortality)
  
}


odesys_mean <- function (t, rho, parms) {
  list( d_rho_mean(rho, parms )  )
}


## pair-approximation


d_rho_1 <- function(rho, parms) { 
  
  growth <- with(parms, 
                 r * (b + (1 - b) * f * (rho[1] - rho[2])/(1- rho[1]) ) * rho[1]^( 1 + alpha) * (1 - rho[1]/(K * (1-c*(rho[1] - rho[2])/(1- rho[1])) ) )) 
  if(growth <= 0 | is.na(growth)) {growth <- 0}

  mortality <- with(parms, m * rho[1] + ( (a + v*rho[2]/rho[1]) * rho[1]^( 1 + q) * L * (1 - p * rho[2]/rho[1]))/(1 +(a + v*rho[2]/rho[1]) * h  * rho[1]^( 1 + q)) )
  if(mortality <= 0 | is.na(mortality)) {mortality <- 0}
  
  return(growth - mortality)
}
  

d_rho_11 <- function(rho,  parms) { 
  
  growth <- with(parms,
       2* (rho[1] - rho[2]) * r * (b + (1 - b) * f * (rho[1] - rho[2])/(1- rho[1]) ) * rho[1]^( 1 + alpha) * (1 - rho[1]/(K * (1-c*(rho[1] - rho[2])/(1- rho[1])) ) ) / (1-rho[1]))
  if(growth <= 0 | is.na(growth)) growth <- 0
  
  mortality <- with(parms, 2 * rho[2] * m  + 2 * rho[2] * ( (a + v*rho[2]/rho[1]) * rho[1]^( 1 + q) * L * (1 - p * rho[2]/rho[1]))/(1 +(a + v*rho[2]/rho[1]) * h  * rho[1]^( 1 + q))  )
  if(mortality <= 0 | is.na(mortality)) mortality <- 0
  
  return(growth - mortality)
}


odesys_spex <- function(t, rho, parms = model_parms) { 
  if(rho[1] < 1e-7 ) {
    rho <- ini_rho(0,0)
  } 
  
    rho_1 <- d_rho_1(rho, parms)
    rho_11 <- d_rho_11(rho, parms)
 
  list(c(
    rho_1 = rho_1,
    rho_11 = rho_11,
    #changes of other pairs and singletons are not calculated
    rho_10 = 0, 
    rho_00 = 0,
    rho_0 = 0
  )  )
}

# function to get initial rho values for pair-approximation

ini_rho <- function(rho_1, rho_11 = NULL) {
  if(is.null(rho_11[1])) out <- {
    c(
      rho_1 = rho_1,
      rho_11 = rho_1*rho_1,
      rho_10 = rho_1*(1-rho_1)*2,
      rho_00 = (1-rho_1)*(1-rho_1),
      rho_0 = 1-rho_1
    ) 
  } else {
    out <- c(
      rho_1 = rho_1,
      rho_11 = rho_11,
      rho_10 = rho_1-rho_11,
      rho_00 = 1-2*rho_1+rho_11,
      rho_0 = 1-rho_1
    )
    if(any(out < 0)) out <- c(
      rho_1 = NA,
      rho_11 = NA,
      rho_10 = NA,
      rho_00 = NA,
      rho_0 = NA
    )
  }
  
  return(out)
}


# function for local growth
# takes global cover, rho_1 and local cover, q_01, as parameters, plus a list object with model parameters. 

G <- function(rho_1, q_01 = "auto", parms = defparms, set = list(NA)) {
  
  if(q_01[1] == "auto") q_01 <- 1
  
  
  parms.names <- names(parms)
  set.names <- names(set)
  m.names <- sort(unique(c(parms.names, set.names)))
  
  parms_temp <- sapply(m.names, function(i) {
    if (i %in% set.names) set[[i]]
    else parms[[i]]
  }, simplify = FALSE)
  
  out <- with(parms_temp, r*(b + (1-b)*f*q_01)*rho_1^(1+alpha)*(1-rho_1/(K * (1-c * q_01) ) ))
  
  out[out < 0] <- 0
  return(out)
} 


# function for local grazing mortality
# takes global cover, rho_1 and local cover, q_11, as parameters, plus a list object with model parameters. 

C <- function(rho_1, q_11 = 1, parms = defparms, set = list(NA)) {
  
  parms.names <- names(parms)
  set.names <- names(set)
  m.names <- sort(unique(c(parms.names, set.names)))
  
  parms_temp <- sapply(m.names, function(i) {
    if (i %in% set.names) set[[i]]
    else parms[[i]]
  }, simplify = FALSE)
  
  with(parms_temp, (m * rho_1) +( (a+v*q_11  ) *rho_1^(1+q)*L*(1-p*q_11))/(1+(a+v*q_11)*h*(rho_1^(1+q)) ) )

}


colpal <- list(mort = c("#000000","#00000040","#00000030","#00000020"), grow = c("#009933","#00993340","#00993330","#00993320") )


# default plot for graphs over G or C

defplot <- function(...) plot(..., 
                              type = "l", 
                              bty = "l", cex = 0.7, las = 1,
                              xlim = c(0,1), ylim = c(0,0.25), 
                              xaxs = "i", yaxs = "i", 
                              xaxp = c(0,1,2), yaxp = c(0,0.25,2), 
                              xlab = "vegetation cover") 



runODE_spex <- function(starting, model_parms, times = c(0,1000))  {
  
    out <- as.data.frame(ode(starting, func = odesys_spex, times = times, parms = model_parms))
    
    names(out) <- c("time", "rho_1", "rho_11", "rho_10", "rho_00", "rho_0")
    
    out$rho_10 <- out$rho_1- out$rho_11
    out$rho_00 <- 1-2*out$rho_1+out$rho_11
    out$rho_0  <- 1- out$rho_1
    
  return(round(out,6))
}

# graphical visualisation of attractor
# simulates trajectories for many differently clustered starting conditions. 
# running pairapprox = TRUE requires parallel backend

attractor <- function(model_parms, rho_1_ini = seq(0,1, length = 41), rho_11_ini = seq(0,1, length = 11), meanfield = TRUE, pairapprox = TRUE, localvals = FALSE) {
  
  defplot(NA,NA, ylab = "plant mortality/growth")
  rho <- seq(0,1,length = 100)
  
  if(meanfield) {
    
    lines(rho,C(rho, rho, model_parms), col = colpal$mort[1], lwd = 2)
    lines(rho,G(rho, rho, model_parms), col = colpal$grow[1], lwd = 2)
    
    
    runmodel_high <- as.data.frame(ode(y = 0.99, func = odesys_mean, times = c(1,1000), parms = model_parms))
    
    points(runmodel_high[2,2],G(runmodel_high[2,2], runmodel_high[2,2], model_parms), xpd = TRUE, pch = 20, cex = 2)
    
    runmodel_low <- as.data.frame(ode(y = 0.0001, func = odesys_mean, times = c(1,1000), parms = model_parms))
    
    points(runmodel_low[2,2],G(runmodel_low[2,2],runmodel_low[2,2],model_parms), xpd = TRUE, pch = 20, cex = 2)
    
    lo <- runmodel_low[2,2]
    hi <- runmodel_high[2,2]
    
    for(i in 1:10) {
      
      mid <- (lo+hi)/2 
      runmodel_mid <- as.data.frame(ode(mid, func = odesys_mean, times = c(1,1.2), parms = model_parms))
      
      if(runmodel_mid[2,2] > mid) { hi <- mid } else { lo <- mid}
      
    }
    
    
    points(mid,G(mid,mid,model_parms), xpd = TRUE, pch = 21, cex = 1.5, bg = "white")
        
    
  }
  
  if(pairapprox) {
    
    ini <- list(
      rho_1 = rho_1_ini,
      rho_11 = rho_11_ini,
      rho_10 = NA,
      rho_00 = NA,
      rho_0 = NA
    )
    
    ini <- expand.grid(ini)
    
    for(x in 1:nrow(ini)) {
      temp <- ini_rho(ini[x,]$rho_1, ini[x,]$rho_11)
      ini[x,]$rho_1 <- temp[1]
      ini[x,]$rho_11 <- temp[2]
      ini[x,]$rho_10 <- temp[3]
      ini[x,]$rho_00 <- temp[4]
      ini[x,]$rho_0 <- temp[5]
    } 
    
    
    ini$m_ini <- C(ini$rho_1, ini$rho_11/ini$rho_1, model_parms)
    ini$g_ini <- G(ini$rho_1, (ini$rho_1-ini$rho_11)/(1-ini$rho_1), model_parms)
    
    ini <- subset(ini, !is.na(m_ini) & !is.na(ini$g_ini) &  ini$g_ini >= 0)
    ini <- cbind(ID = 1:nrow(ini),ini)
    
    
    foreach(iteration = ini$ID, .combine = rbind, .packages = c("deSolve")) %dopar% { 
      source("C:/Users/SCHNEIDER/Documents/projects/CAS02_livestock/code/simfunctions.r")
      
      rho_starting <- ini[iteration, 2:6]
      
      # running the ode-solver
      runmodel <- runODE_spex(as.numeric(rho_starting),  model_parms, times = seq(1,2))
      
      return(tail(runmodel,1))
    } -> output 
    
    steady <- ini$rho_1 == output$rho_1
    ini <- ini[!steady,]
    output <- output[!steady,]
    
    arrows(ini$rho_1, C(ini$rho_1, ini$rho_11/ini$rho_1, model_parms),output$rho_1, C(output$rho_1, output$rho_11/output$rho_1, model_parms), length = 0.02)
    
    arrows(ini$rho_1, G(ini$rho_1, (ini$rho_1-ini$rho_11)/(1-ini$rho_1), model_parms),output$rho_1, G(output$rho_1, (output$rho_1-output$rho_11)/(1-output$rho_1) , model_parms), length = 0.02, col = "#009933" )
    
    ## plot steady states
    
    high_equ <- as.data.frame(ode(y = ini_rho(0.8), func = odesys_spex, times = c(1,1000), parms = model_parms) )[2,]
    
    points(high_equ$rho_1, C(high_equ$rho_1, high_equ$rho_11/high_equ$rho_1, model_parms), pch = 20, xpd = TRUE)
    
    low_equ <- as.data.frame(ode(y = ini_rho(0.0001), func = odesys_spex, times = c(1,1000), parms = model_parms) )[2,]
    
    points(low_equ$rho_1, C(low_equ$rho_1,low_equ$rho_11/low_equ$rho_1, model_parms), pch = 20, xpd = TRUE)
    
  }
  
  if(localvals) {
    
 
    lines(rho,C(rho, 1, model_parms), col= colpal$mort[3])
    lines(rho,C(rho, 0.75, model_parms), col= colpal$mort[3])
    lines(rho,C(rho, 0.5, model_parms), col= colpal$mort[2])
    lines(rho,C(rho, 0.25, model_parms), col= colpal$mort[3])
    lines(rho,C(rho, 0, model_parms), col= colpal$mort[4])
    
    lines(rho, G(rho, 1, model_parms), col = colpal$grow[4])
    lines(rho,G(rho, 0.75, model_parms), col = colpal$grow[3])
    lines(rho,G(rho, 0.5, model_parms), col = colpal$grow[2])
    lines(rho,G(rho, 0.25, model_parms), col = colpal$grow[3])
    lines(rho,G(rho, 0, model_parms), col = colpal$grow[4])
    
  }
  
  
}


bifurcation <- function(parms, over, xrange, res = 201, times = c(0,1000), ini = c(0.9, 0.0001) , meanfield = TRUE, pairapprox = FALSE, numerical = FALSE ) {
  
  require(foreach)
  
  parms[[over]] <- seq(xrange[1],xrange[2],length = res)
  
  parms$rho_ini <- ini
  
  iterations <- expand.grid(parms)
  iterations <- cbind(ID = 1:dim(iterations)[1],iterations)
      
  foreach(iteration = iterations$ID, .combine = rbind, .packages = c("deSolve")) %dopar% { 
    source("C:/Users/SCHNEIDER/Documents/projects/CAS02_livestock/code/simfunctions.r")
    
    model_parms <- as.list(iterations[iteration,])
    
    # running the ode-solver
    runmodel <- runODE_spex(ini_rho(model_parms$rho_ini), model_parms, times = times) 
      
    return(tail(runmodel,1))
  } -> output
  
  output <- cbind(iterations,output)
  
  upper <- output[output$rho_ini == ini[1],][which(round(output[output$rho_ini == ini[1],]$rho_1,4) != round(output[output$rho_ini == ini[2],]$rho_1,4)),]
  lower <- output[output$rho_ini == ini[2],][which(round(output[output$rho_ini == ini[2],]$rho_1,4) != round(output[output$rho_ini == ini[1],]$rho_1,4)),]
  
  if(nrow(upper)>0) {
  foreach(i = upper[,over], .combine = rbind, .packages = c("deSolve") ) %dopar% {
    
    source("C:/Users/SCHNEIDER/Documents/projects/CAS02_livestock/code/simfunctions.r")
    
    model_parms <- upper[upper[, over] == i,]
    
    hi_1 <- upper[upper[, over] == i,]$rho_1
    lo_1 <- lower[lower[, over] == i,]$rho_1
    hi_11 <- upper[upper[, over] == i,]$rho_11
    lo_11 <- lower[lower[, over] == i,]$rho_11
    
    for(j in 1:10) {
      
      rho_ini <- ini_rho( (hi_1+lo_1)/2 , (hi_11+lo_11)/2 )
      
      # running the ode-solver
      
      runmodel <- runODE_spex(rho_ini, model_parms,times = c(0,1.5)) 
      
      if(runmodel[2,"rho_1"] < runmodel[1,"rho_1"] ) {
        lo_1 <- (hi_1+lo_1)/2 
        lo_11 <- (hi_11+lo_11)/2 
      } else {
        hi_1 <- (hi_1+lo_1)/2 
        hi_11 <- (hi_11+lo_11)/2 
      }
      
    }
    
    return(tail(runmodel,1))
  } -> output_unstable
  
  
  output_unstable <- cbind(upper[,1:16],output_unstable)
  
  }
  
  plot(output$rho_1 ~ output[,over], xlab = over, ylab = "vegetation cover", type = "p", pch  = 20, ylim = c(0,1), cex = 0.5, yaxp = c(0,1,2))
  
  if(nrow(upper)>0) {
  points(output_unstable$rho_1 ~ output_unstable[,over], pch = 20, col = "grey80", cex = 0.5)
}

}
  
  