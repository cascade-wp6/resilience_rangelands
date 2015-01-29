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
source("code/simfunctions.r")

set.seed(584232) # setting seed for random number generation

# parameter settings:
parameters <- list(  
  m = 0.05, #intrinsic mortality
  r = 0.5, # intrinsic growth rate
  b = 1,   # environmental quality / 1 - aridity
  K = 0.9, # carrying capacity 
  a = 0.6, # search efficiency
  h = 100, # handling time
  L = 10 # livestock units / grazing intensity
)


nullmodel <- runCA(runif(1, 0.8,0.9), parameters, delta = 0.2, t_max = 150, t_min = 100, t_eval = 50, stability = 0.0001, saveeach = 5)


par(mfrow = c(2,1))
plot(nullmodel$rho_one[0:i] ~ nullmodel$time[0:i], 
     ylim = c(0,1), ylab = expression( rho["+"]), 
     xlab = "time [yr]", type = "l")
     

plot(nullmodel$q_one_one[0:i] ~ nullmodel$time[0:i], 
     ylim = c(0,1), ylab = expression(bar(q["1|1"] )), 
     xlab = "time [yr]", type = "l")


animateCA(nullmodel, "ca_nullmodel.gif")
