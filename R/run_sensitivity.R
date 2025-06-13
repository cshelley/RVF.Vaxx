### Package installations
if (!require("sensitivity")) {       #Sensitivity Analysis Package
  install.packages("sensitivity", dependencies = TRUE)
  library(sensitivity)
}

if (!require("FME")) {               #Sensitivity Analysis Package
  install.packages("FME", dependencies = TRUE)
  library(FME)
}


### Model Sensitivity Analysis ###


## 1. Generate Latin Hypercube Sample


vars <- rbind(nuH = c(0,1), rhoVH = c(0,1), rhoRH = c(0,1),
              betaAH = c(0,1), alphaH = c(0,1),
              nuM = c(0,1), rhoVM = c(0,1), rhoRM = c(0,1),
              betaAM = c(0,1), alphaM = c(0,1),
              tau = c(0,1), omega = c(0,1), zeta = c(0,1),
              betaHA = c(0,1), betaMA = c(0,1), muA = c(0,1), muE = c(0,1))


### This is the hypercube of all my multipliers with [0,1] ranges
### (the extremes of what is possible within my model)

factors <- c("nuH", "rhoVH", "rhoRH", "betaAH", "alphaH", "betaAM", "alphaM",
             "tau", "omega", "zeta", "betaHA", "betaMA", "muA", "muE")


parRange <- data.frame(min = vars[,1], max = vars[,2])
#rownames(parRange) <- factors      ### Type "parRange" to see the resulting data.frame


LHS <- data.frame(Latinhyper(parRange, 1000))
## Latin Hypercube sampling, 1000 possible combinations of parameter values


LHSplot <- pairs(LHS, main = "Latin Hypercube")
## A picture of the LHS, a variance-covariance matrix. Try rerunning LHS on line 34 to be
## only 10 resamplings and you'll see well what this is doing. It's dividing each covariance
## plot into a 10x10 grid and choosing a point within each column and each row. It leaves serious
## blank spots which is why Sobol' sampling is better but I could never get that to run properly.


init = c(VH = 0, SH = 10000, IH = 2, RH = 0,
         SA = 175000, IA = 100, QA = 0, PA = 0,
         VM = 0, SM = 100000, IM = 200, RM = 0)

times = 1:500

y <- vector()
for(i in 1:nrow(LHS)) {
  rvf <- rvf_sensitivity(t = times, state = init, parameters = LHS[i,])
  y <- append(y, rvf, after = length(y))
}


## Here I am running my model using each row of the LHS as my parameter values.
## My model currently takes parms as a vector of inputs so I'm just switching to
## 1000 runs of LHS[i,] in lieu of parms


PCC <- pcc(X = LHS, y = unlist(y), rank = TRUE, nboot = 100, conf = 0.95)
## This is the PRCC calculation (seriously, that's it). My confidene is set to 0.95
## For a simple calculation you can just put in X and it will do the math itself, but
## my model is doing too much for the command to handle and so I'm also providing it
## with y = outcome to x = income. nboot is the number of bootstrap resampling's you're
## doing. You need this to come up with the confidence limits. Don't make it too high or
## it'll run all day. I tried a few values and saw that 100 is stable but small enough to run.


## Type "PCC" to see what you actually generated here. That's the order of parameters
## you'll see in the plot coming up, but the auto-plot won't give x-axis labels so I'm
## doing them by hand.




## Plot


prcc <- PCC$PRCC[,1]
print(PCC$PRCC)
barplot(prcc)   #x axes are in the same order as print(PRCC$PRCC) above


## Publication quality plot, make sure x axes match pcc$PRCC print out
par(mfrow = c(1,1), cex = 1.2)
prcc <- PCC$PRCC[,1]
barplot(prcc, ylim = c(-0.2, 0.2), xaxt='n', ann=FALSE, yaxt='n',
        xlab = "Parameter Value", ylab = "PRCC")
abline(h=0, col = "red")
xes <- c(0.7, 1.9, 3.1, 4.3)
segments(x0=xes, x1=xes, y0 = PCC$PRCC[,4], y1=PCC$PRCC[,5])        # vertical whiskers
segments(x0=xes-.1, x1=xes+.1, y0 = PCC$PRCC[,4], y1=PCC$PRCC[,4])  # bottom horizontal
segments(x0=xes-.1, x1=xes+.1, y0 = PCC$PRCC[,5], y1=PCC$PRCC[,5])  # top horizontal


axis(1, at = xes, labels = expression(beta, epsilon, gamma, kappa))


axis(2, at = seq(-0.2, 0.2, by = 0.05), label = seq(-0.2, 0.2, by = 0.05))
box()



################################
# Cut:

## Deterministic measles function
measles <- function(times, init, parms) {
  derivs <- function(t, init, parms) {
    with(as.list(c(init, parms)), {


      sigma.b <- kappa*sqrt(beta*S*exp(log(I)))/sqrt(N)
      sigma.e <- kappa*sqrt(epsilon*exp(log(E)))/sqrt(N)
      sigma.g <- kappa*sqrt(gamma*exp(log(I)))/sqrt(N)

      # stochastic noise terms
      zeta.b = rnorm(1,0, sigma.b)
      zeta.e = rnorm(1,0, sigma.e)
      zeta.g = rnorm(1,0, sigma.g)

      # SDE system of equations
      #dM <- -theta1*M
      dS <- -beta*exp(log(S))*exp(log(I)) - zeta.b
      dE <- beta - epsilon*E  + zeta.b - zeta.e
      dI <- epsilon*E - gamma + zeta.e - zeta.g
      dR <- gamma*exp(log(I)) + zeta.g
      #dV1 <- omega1*S-theta2*V1-omega2*V2
      #dV2 <- omega2*V2-theta3*V2

      #list(c(dM,dS,dE,dI,dR,dV1,dV2))
      list(c(dS,dE,dI,dR))
    })
  }
  result <- ode(y=init, times = times, func = derivs, parms=parms)
  as.data.frame(result)
}



