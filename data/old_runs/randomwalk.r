x_new 


grazerlayer <- function(x, n, state = c("0", "-"), start = "random") {

starting <- sample( which(x$cells %in% state & count(x_new,  c("0", "-")) > 3 ), n) 
Y <- matrix(FALSE, ncol = x$dim[1], nrow = x$dim[2]) 
t(Y)[starting] <- TRUE
out <- list()
out$pos <- starting
out$arr.ind <- which(Y, arr.ind = TRUE)
out$cumulative <- Y
out$orientation <- numeric(length = n)
out$saturation <- numeric(length = n)
out$state <- state
class(out) <- c("list", "grazerlist")
return(out)
}

out <- grazerlayer(x_new, 50)

par(mar =c(0,0,0,0))
plot(x_new, cols = color)
points(out$arr.ind, pch = 21, bg = "red3")

move <- function(x, steps, landscape) {
n = length(x$pos)
opencells_logical_with_border <- (landscape$cells %in% x$state)[x_with_border]

for(i in 1:steps) {
# for each x$pos get  a vector of type allow <- c(TRUE, TRUE, FALSE, TRUE ) for allowance to walk east,  north, west, south
allow <- sapply(interact, function(k) opencells_logical_with_border[x_to_evaluate[x$pos]+k])

# for each x$pos sample a random value from  c(1,2,3,4) with probability allow 
goto <- sapply(1:n , function(k) sample(c(1:4), 1, prob = allow[k,]) )

# calculate new x$pos
newind <- x$arr.ind + data.frame(row = c(-1,0,1,0)[goto], col = c(0,1,0,-1)[goto])
newind$row[newind$row > landscape$dim[1]] <- newind$row[newind$row > landscape$dim[1]]- landscape$dim[1]
newind$col[newind$col > landscape$dim[2]] <- newind$col[newind$col > landscape$dim[2]]- landscape$dim[2]

Y <- matrix(FALSE, ncol = landscape$dim[1], nrow = landscape$dim[2]) 

set <- function(x,y) Y[x, y] <- TRUE
set(newind$row, newind$col)

which(Y == TRUE)

x$cumulative 

points(grazerlist, pch = 21, bg = "red3")
}

}



