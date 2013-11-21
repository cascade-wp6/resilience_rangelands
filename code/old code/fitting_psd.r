
setwd("E:\\Eigene Dokumente\\Uni\\projects\\2013 CASCADE grazing\\sim4")

output <- read.csv("output.csv")

output[(output$b == 0.52 & output$globalgrazing == TRUE & output$starting == "high"),]

load(paste("results\\result", 777, sep = "_"))





pdf("E:\\SkyDrive\\Uni\\talks\\patchdist_graz.pdf", height = 5, width = 6, paper = "special")
par(mfrow = c(1,1))

for(j in seq(0,0.025, 0.005)) {
	plot(NA, NA, xlim = c(1,10000), ylim = c(1,1000), log= "xy", main = paste("g =", j))

    select <- with(output, ID[which(g == j & b == 0.52 & globalgrazing == TRUE & starting == "high")])
	load(paste("results/result_", select, sep = ""), envir = globalenv() )
	if(all(!is.na(result$cumpatch))) points(mean~size, pch = 20, col = "#968376", data = result$cumpatch)

	
    select <- with(output, ID[which(g == j & b == 0.52 & globalgrazing == FALSE & starting == "high")])
	load(paste("results/result_", select, sep = ""), envir = globalenv() )
	if(all(!is.na(result$cumpatch))) points(mean~size, pch = 20, col = "blue3", data = result$cumpatch)


}
dev.off()


















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




with(result$cumpatch, 
	plot(log(mean)~log(size), type = "p", pch = 20)
	
	dd1 <- data.frame(size = NA, n = NA )[-1,]
	for(i in as.character((snapshots-1)*delta)[-(1:2)] ) {
	dd1 <- rbind(dd1, data.frame(size = result$cumpatch[,1], n = result$cumpatch[,i] ))
	}
	
	
)

plot(mean~size, type = "p", pch = 20, data = result$cumpatch, log = "xy")


model1 <- nls(n ~ a * size^I(-alpha), 
		data = dd1,
		start = list(a = 47.760, alpha = 1.5000),
        trace = TRUE
		)
		
lines(	exp(seq(log(1), log(200), length = 100)), 
		predict(model1, list(size = exp(seq(log(1), log(200), length = 100   )))) 
	)
	

model2 <- nls(n ~ a * size^I(-alpha) * exp(-(size/Sx)), 
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
	