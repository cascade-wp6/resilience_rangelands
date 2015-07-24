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

rm(list=ls())

########################################################################################
library(deSolve)
source("code/simfunctions.r")

set.seed(584232) # setting seed for random number generation


# default parameters for spatially explicit model

defparms <- list(
  m = 0.05, # intrinsic mortality
  r = 1,    # maximal growth rate
  b = 0.7,  # environmental quality
  K = 0.9,  # carrying capacity
  a = 0.2,  # attack rate of livestock
  h = 50,   # handling time of livestock in [time per area]
  L = 4,    # number of livestock
  alpha = 0,# water runoff
  q = 0,    # hill-coefficient of functional response
  c = 0,    # local competition
  f = 0,    # local facilitation
  v = 0,    # local attraction
  p = 0     # local protection
)


# do a single simulation run


## set your parameters
parms <- defparms
parms$c = 1
parms$L = 6
parms$p = 1


# run spatially explicit model

out_spex <- runODE_spex( ini_rho(0.8), model_parms = parms,  times = exp(seq(0,4,length = 100))-1)

plot(rho_1 ~ time , data = out_spex , type = "l", ylim = c(0,1), xlim = c(0,50))
axis(4, at = round(tail(ODE_spex,1)$rho_1,2),cex.axis = 0.7, las = 1 )


# run mean-field model

out_meanfield <- runODE_meanfield( ini_rho(0.8), model_parms = parms,  times = exp(seq(0,4,length = 100))-1)

lines(rho_1 ~ time , data = out_meanfield, lty = 2)


# run cellular automata model

out_ca <- runCA(0.8, parms, delta = 0.1, t_max = 150, t_min = 100, t_eval = 50, isstable = 0.0001, saveeach = 5)

lines(out_ca$rho_one ~ out_ca$time)

# provide parallel backend

library(foreach)
library(doSNOW)


# this is to run the cluster on your local computer:
workerlist <- c(rep("localhost", times = 7)) 

cl <- makeSOCKcluster(workerlist, outfile='out_messages.txt')

# THIS IS NOT GOING TO WORK ON THE WORKSTATION --------------------
# if you need to use the workstation specify the following 
#workstation <-  list(host = "162.38.184.118",  # IP of the host
#                     user = "florian",  #change to your username
#                     rscript = "/usr/lib/R/bin/Rscript", # should be correct, i.e. the place where  R is installed on the workstation 
#                     snowlib = "/usr/lib/R/library")  # the place where SNOW package is installed on the workstation, you might need to install doSNOW on the workstation with admin rights (login to your account, start 'sudo R', then 'install.packages("doSNOW")')

#workerlist <- c(rep(workstation, times = 23)) 

#cl <- makeSOCKcluster(workerlist, master="162.38.184.xx", outfile='out_messages.txt') # add your computer's IP address

# THIS IS NOT GOING TO WORK ON THE WORKSTATION --------------------

registerDoSNOW(cl)



# visualize attractor (isoclines)

## meanfield 

attractor(parms, pairapprox = FALSE, localvals = TRUE)


### .----- not working beyond this point --- code is too specific for my computer --- sorry --- need to wrap it into a package ---


## pair-approximation

attractor(parms, pairapprox = TRUE, meanfield = FALSE, localvals = FALSE)


# create bifurcation diagram

## meanfield over b
bifurcation(parms, pairapprox = FALSE, "b", c(0,1))

## meanfield over L

bifurcation(parms, pairapprox = FALSE, "L", c(0,12), res = 51)

## pair-approximation 
bifurcation(parms, pairapprox = TRUE, meanfield = FALSE, "b", c(0,1))






stopCluster(cl)

