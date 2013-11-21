size <- exp(seq(log(1), log(10000), length = 200))

a <- 0.8
alpha <- 1
b <- .01
Sx <- 100

PL <-  exp(  log(a) - alpha * log(size) )
TPLup <-  exp(  log(a) - alpha * log(size) + log(1+b/(a*size^(-alpha)))) 
TPLdown <-  exp(    log(a) - alpha * log(size) - (size/Sx)  )

eps = .3

EXP <-   exp( log(a) - (eps*size) )
plot(NA, NA, xlim = c(1, 10000) , ylim = c(0.0001, 1), log = "xy")
lines(size, PL, lwd = 2)
lines(size, TPLup)
lines(size, TPLdown)
lines(size, EXP)


