
# provides parallel backend
library(foreach)
library(doSNOW)

source("code/simfunctions.r")
workstation <-  list(host = "162.38.184.118", user = "schneider",
                     rscript = "/usr/lib/R/bin/Rscript",
                     snowlib = "/usr/lib/R/library")

workerlist <- c(rep("localhost", times = 7)) 

cl <- makeSOCKcluster(workerlist, outfile='out_messages.txt')

registerDoSNOW(cl)


library("deSolve")

defaultparms <- list(
  m = 0.05,  # intrinsic mortality of plants 
  r = 1,   # growth rate
  f = 0, #  local facilitation
  b = 0.85,  # environmental quality inverse of aridity
  K = 0.9,  # carrying capacity
  alpha = 0, # water runoff
  c = 0, # local competition
  a = 0.3, # search efficiency
  v = 0, # attraction effect
  h = 50, # handling time (one individual on landscape unit)
  p = 0, # protection effect
  L = 2, # livestock units per landscape unit 
  q = 0 # search ineffiency at low cover
)




parms <- defaultparms
parms$m = 0
parms$L = 0
parms$b = 0
parms$f = 1
l0 <- runCA(0.005, parms, 25,25, t_max = 400)



parms <- defaultparms

l1 <- runCA(0.05, parms, 25,25)
parms$L <- 8
l2 <- runCA(0.95, parms, 25,25)

png("talk_fig_2.png", width = 4, height = 6, units = "in", res = 300)

par(mfrow = c(3,2), mar = c(1,1,1,1))
plot(l1$timeseries[[1]])
box()
plot(l1$timeseries[[8]])
box()
plot(l2$timeseries[[1]])
box()
plot(l2$timeseries[[4]])
box()
plot(l0$timeseries[[1]])
box()
plot(l0$timeseries[[6]])
box()
dev.off()



parms <- defaultparms

png("talk_fig_1a.png", width = 8, height = 3.5, units = "in", res = 300)
par(mfrow = c(1,2), bty = "n", las = 1)
attractor(parms, localvals = TRUE, pairapprox = FALSE, meanfield = TRUE, rho_1_ini = seq(0,1, length = 11) )

bifurcation(parms,"L", c(0,12)) -> over_b
abline(v = 2)
dev.off()

png("talk_fig_1b.png", width = 8, height = 3.5, units = "in", res = 300)
parms$L = 4
par(mfrow = c(1,2), bty = "n", las = 1)
attractor(parms, localvals = TRUE, pairapprox = FALSE, meanfield = TRUE, rho_1_ini = seq(0,1, length = 11) )

bifurcation(parms,"L", c(0,12)) -> over_b
abline(v = 4)
dev.off()

png("talk_fig_1c.png", width = 8, height = 3.5, units = "in", res = 300)
parms$L = 8
par(mfrow = c(1,2), bty = "n", las = 1)
attractor(parms, localvals = TRUE, pairapprox = FALSE, meanfield = TRUE, rho_1_ini = seq(0,1, length = 11) )

bifurcation(parms,"L", c(0,12)) -> over_b
abline(v = 8)
dev.off()




png("talk_fig_1d.png", width = 8, height = 3.5, units = "in", res = 300)
parms$L = 8
parms$p = 1
parms$a = 1
par(mfrow = c(1,2), bty = "n", las = 1)
attractor(parms, localvals = TRUE, pairapprox = FALSE, meanfield = TRUE, rho_1_ini = seq(0,1, length = 11) )

bifurcation(parms,"L", c(0,12)) -> over_b
abline(v = 8)
dev.off()


png("talk_fig_1e.png", width = 8, height = 3.5, units = "in", res = 300)
parms$L = 8
parms$p = 1
parms$a = 1
par(mfrow = c(1,2), bty = "n", las = 1)
attractor(parms, localvals = TRUE, pairapprox = TRUE, meanfield = FALSE, rho_1_ini = seq(0,1, length = 11) )

bifurcation(parms,"L", c(0,12), pairapprox = TRUE) -> over_b
abline(v = 8)
dev.off()
