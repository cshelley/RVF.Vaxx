library(sensitivity)
library(FME)
library(parallel)

## 1. Generate Latin Hypercube Sample
vars = rbind(nuH = c(0,1), omega = c(0,1), alphaH = c(0,1), nuM = c(0,1),
             muA = c(0,1), alphaM = c(0,1), tau = c(0,1), rhoRM = c(0,1),
             betaM = c(0,1), rhoRH = c(0,1), rhoVM = c(0,1), muE = c(0,1),
             betaAM = c(0,1), betaMA = c(0,1), betaHA = c(0,1), rhoVH = c(0,1),
             zeta = c(0,1), betaAH = c(0,1))

parRange <- data.frame(min = vars[,1],
                       max = vars[,2])
### This is the hypercube of all my multipliers with [0,1] ranges
### (the extremes of what is possible within my model)


LHS <- data.frame(Latinhyper(parRange, 5000))
## Latin Hypercube sampling, 10000 possible combinations of parameter values

#LHSplot <- pairs(LHS, main = "Latin Hypercube")
## A picture of the LHS, a variance-covariance matrix. Try rerunning LHS to be
## only 10 resamplings and you'll see well what this is doing. It's dividing
## each covariance plot into a 10x10 grid and choosing a point within each
## column and each row. It leaves serious blank spots which is why Sobol'
## sampling is better but I could never get that to run properly.


## 2. Run LHS parameters through disease model to generate
init = c(SH = 1000, IH = 0, RH = 0, VH = 100,
          SM = 10000, IM = 20, RM = 0, VM = 100,
          SA = 100000, IA = 200, PA = 0, QA = 0)


times = 1:50

y <- vector()
for(i in 1:nrow(LHS)) {
  rvf <- ode(y = init, times = times, func = rvf_sensitivity, parms = LHS[i,])
  y <- append(y, rvf[50,3], after = length(y))
}

## Here I am running my model using each row of the LHS as my parameter values.
## My model currently takes parms as a vector of inputs so I'm just switching to
## 1000 runs of LHS[i,] in lieu of parms

## I'm also only really interested in human infections so I'm only saving number
## of infections at the 50th time step

## This is pretty slow so I should parallelize this step.


PCC <- pcc(X = LHS, y = unlist(y), rank = TRUE, nboot = 1000, conf = 0.95)
## This is the PRCC calculation (seriously, that's it).
## My confidence is set to 0.95.
## For a simple calculation you can just put in X and it will do the math
## itself, but my model is doing too much for the command to handle and so I'm also providing it
## with y = outcome to x = income. nboot is the number of bootstrap resampling's you're
## doing. You need this to come up with the confidence limits. Don't make it too high or
## it'll run all day. I tried a few values and saw that 100 is stable but small enough to run.


## Type "PCC" to see what you actually generated here. That's the order of parameters
## you'll see in the plot coming up, but the auto-plot won't give x-axis labels so I'm
## doing them by hand.
