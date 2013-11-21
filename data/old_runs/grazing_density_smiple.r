
# specify lattice
width = 100
height = 100

states = c("+","0","-")
color <- c("black","grey80", "white") # define colors for the cell state levels

parms_temp <- list(
	m = 0.1,
	g = .1

)


res = 200
death <- data.frame(dens = seq(0,1, length=res))

for(i in 1:res) {

j = i/res 
prob = c( i/res ,1- i/res ,0)

#### first run with grazing
# sampling the initial random grid into a list object
initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = prob ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)


parms_temp$rho  <-  sum(initial$cells == "+")/(width*height)
parms_temp$Q_plus <- count(initial, "+")/4   
parms_temp$vul <- sum((1-parms_temp$Q_plus)[initial$cells == "+"])/(width*height)


#if(parms_temp$vul == 0) parms_temp$vul = 1/(width*height)

death$parent[i] <- with(parms_temp, (m+g/0.75))
death$stock[i] <- with(parms_temp, (m+g/rho))
death$assoc[i] <- mean(with(parms_temp, (m+g*(1-Q_plus)/0.1888)))
death$stockassoc[i] <- mean(with(parms_temp, (m+g*(1-Q_plus)/vul)))


}

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), ylab = "risk of death", xlab = "vegetation cover")

lines(parent*dens ~ dens, data = death, col = "black", lwd = 2)

lines(stock*dens ~ dens, data = death, col = "red3", lwd = 3)
lines(assoc*dens ~ dens, data = death, col = "blue", lwd = 2)
lines(stockassoc*dens ~ dens, data = death, col = "green2", lwd = 1)
