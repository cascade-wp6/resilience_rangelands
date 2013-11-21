
#setwd("/home/schneider/herbivory/")
setwd("C:\\Users\\SCHNEIDER\\Documents\\CASCADE\\2013_grazingmodel\\corrected_grazing_vul\\sim4")

load(paste("results\\result", 638, sep = "_"))

sort(unlist(result$patches))

result$cumpatch$size
plot(mean~size, data = result$cumpatch, log= "xy", pch = 20)

lmean <- mean(log(sort(unlist(result$patches))))
lstd  <- sd(log(sort(unlist(result$patches))))

logncdf <- function(x, lmean, lstd) 1/(std*sqrt(2*pi)*sum(exp(-(log())) seq())) 
P <- ecdf(sort(unlist(result$patches)))

qlnorm(sort(unlist(result$patches)))



	bins <- sort(unique(unlist(result$patches[-1])))
	bins <- seq(1,501, 10)
	bins <- 10^seq(log10(1),log10(10000), .1)
	plot(NA,NA, xlim = c(1,10000), ylim = c(1,1000), log = "xy")
	for(i in 3:22) {
	freq <- sapply(bins, function(j) length(which(result$patches[[i]] >= j)) )
	
	points(freq~bins, pch = 20, type = "b")
}

	with(result$cumpatch2, lines(mean~size, lwd = 3))
	with(result$cumpatch2, lines(c(I(mean+sd),rev(I(mean-sd)))~c(size, rev(size)), lty = 2))

	bins <- sort(unique(unlist(result$patches[-1])))
	plot(NA,NA, log = "xy", xlim = c(1,1000), ylim = c(1,1000))
	for(i in 3:22) {
	freq <- sapply(bins, function(j) length(which(result$patches[[i]] == j)) )
	
	points(freq~bins, pch = 20)
}



# time and resolution of simulation
timesteps = 3000
addgrazing = 1000
delta = 1/5

snapshots <- c(1, c(addgrazing, timesteps-20:1*50, timesteps)/delta+1)

plot(mean~size, type = "p", pch = 20, data = result$cumpatch, log = "xy")


dd1 <- data.frame(size = NA, n = NA )[-1,]

	for(i in as.character((snapshots-1)*delta)[-(1:2)] ) {
	dd1 <- rbind(dd1, data.frame(size = result$cumpatch[,1], n = result$cumpatch[,i] ))
	}
	
ddmean <- data.frame(size = result$cumpatch[,1], n = result$cumpatch$mean )	
	
points(dd1$size, dd1$n, pch = 20, col = "#00000030")

model1 <- nls(n ~ a * size^I(-alpha), 
		data = dd1,
		start = list(a = 47.760, alpha = 1.5000),
        trace = TRUE
		)
		
lines(	exp(seq(log(1), log(200), length = 100)), 
		predict(model1, list(size = exp(seq(log(1), log(200), length = 100   )))) 
	)


model2 <- nls(n ~ a * size^I(-alpha) * exp((size/Sx)), 
		data = dd1,
		start = list(a = 47.760, alpha = 1.5000, Sx = 3),
        trace = TRUE
		)
		
lines(	exp(seq(log(1), log(200), length = 100)), 
		predict(model2, list(size = exp(seq(log(1), log(200), length = 100   )))) 
	)
	
model3 <- nls(n ~ a *  exp(-eps*size) , 
		data = dd1,
		start = list(a = 107.6,eps = 0.8348),
        trace = TRUE
		)
		
lines(	exp(seq(log(1), log(200), length = 100)), 
		predict(model3, list(size = exp(seq(log(1), log(200), length = 100   )))) 
	)		
	
	
	#	model4 <- nls(n ~ 1/(size*sigma*sqrt(2*pi)) * exp((log(size)-mu)^2 / 2*sigma^2)  ,   #from master thesis


model4 <- nls(n ~ 1/2 + (size*sigma*sqrt(2*pi)) * exp((log(size)-mu)^2 / 2*sigma^2)  , 
		data = dd1,
		start = list( mu = 1, sigma = 1),
        trace = TRUE
		)
		
mu = sum(log(dd1$size))/length(dd1$size)
sigma = sqrt(sum((log(dd1$size)- mu)^2) / length(dd1$size))


plot(mean~size, type = "p", pch = 20, data = result$cumpatch, log = "xy", ylim = c(0.1,100))

lines(	exp(seq(log(1), log(200), length = 100)), 
		(1-plnorm( exp(seq(log(1), log(200), length = 100))  ,meanlog = 1-mu, sdlog = sigma))*200
	)
	