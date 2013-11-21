library(rbenchmark)

benchmark(
check1 = {		if(i == addgrazing/delta ) parms_temp <- as.list(parameters[iteration,]) },
replace1 = {		x_new <- x_old 	},
getrho = {		parms_temp$rho <- result$rho[[1]][i-1] },
countQplus = { 		parms_temp$Q_plus <- count(x_old, "+")/4   },
countQempty = {		parms_temp$Q_unveg <- ( count(x_old, c("0", "-")))/4 },
		#parms_temp$vul <- sum(parms_temp$Q_unveg[x_old$cells == "+"]*4)/(width*height)
getrandomnum = {		rnum <- runif(width*height) },
	
calcprobs = {	
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		death <- with(parms_temp, (m+g*Q_unveg)*delta)
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		},
		
checkprobs = {		if(any(c(recolonisation+degradation, death, regeneration) > .95 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
		},
applyrules = {		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"
},

savemort = {
		result$mortality[i] <- mean(death[x_old$cells == "+"], na.rm = TRUE)/delta
		},
		
savebordermort = {		result$mortality_border[i] <- mean(death[x_old$cells == "+" & parms_temp$Q_unveg > 0], na.rm = TRUE)/delta
		},

saverho = {		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}
},

saveqplus2 = {
	result$q_[[1]][i]  <- mean(parms_temp$Q_plus[x_new$cells == "+"], na.rm = TRUE)
},	

savesnapshots = {
	if(i %in% snapshots) write.table(as.data.frame(matrix(c( parms_temp$ID, TRUE, parms_temp$starting, 1, parms_temp$g, parms_temp$m, parms_temp$b, parms_temp$g, parms_temp$m, initial$cells), nrow = 1, dimnames = list("1", header_grids))), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
},

savesnapshots2 = {
if(i %in% snapshots)result$timeseries <- rbind(result$timeseries, 
								as.data.frame(matrix(c( parms_temp$ID, TRUE, parms_temp$starting, 1, parms_temp$g, parms_temp$m, parms_temp$b, parms_temp$g, parms_temp$m, initial$cells), nrow = 1, dimnames = list("1", header_grids)))
								)
},

replace1 = {		x_old <- x_new
		},
		
#clearmem = {		gc()},  #garbage collection
replications = 5000)
