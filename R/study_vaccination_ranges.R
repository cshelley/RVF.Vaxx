# define exploratory ranges for vaccine coverage
domain = 10
human_range = seq(0, 10, length.out = domain)
livestock_range = seq(0, 100, length.out = domain)

# define empty matrices to store output
max_infected = matrix(rep(NA, domain^2), nrow = domain)
time_to_peak = matrix(rep(NA, domain^2), nrow = domain)

# define estimates that won't dynamically define in code
culex_point_estimates = c(betaAH = 0.01, betaAM = 0.01, tauA = 0.2,
                          zetaA = 0.5, omegaA = 0.01, kappaA = 175000,
                          gammaA = 10, omegaA11 = 0.01, muPA = 0.00001,
                          muSA = 1/3, muIA = 1/3, muQA = 0.00001)

aedes_point_estimates = c(betaBH = 0.01, betaBM = 0.01, tauB = 0.2,
                          zetaB = 0.05, omegaB = 0.01, kappaB = 175000,
                          gammaB = 25, omegaB11 = 0.01, muPB = 0.005,
                          muSB = 0.1, muIB = 0.1, muQB = 0.005)

# loop through vaccine ranges and capture time to peak, height of peak
for(nuH in human_range) {
  for(nuM in livestock_range) {

    human_point_estimates =  c(nu_H = nuH/100, nuTH = 0.9, betaHA = 0.81,
                               betaHB = 0.89, epsilonH = 1/4, alphaH = 1/4,
                               rhoRH = 1/900, rhoVH = 1/(6*365), bH = 1/(5*365),
                               muH = 1/(70*365), gammaH = 1.5/(70*360),
                               xiH = 1/4, thetaH = 1/14)

    livestock_point_estimates = c(nu_M = nuM/100, nuTM = 0.9, betaMA = 0.81,
                                  betaMB = 0.89, xiM = 24/3.25,
                                  alphaM = 1/(3*0.95), rhoRM = 1/900,
                                  rhoVM = 1/360, muM = 0.0008, gamma1 = 0.00082,
                                  gamma2 = 0.000082, xiM = 24/3.25,
                                  thetaM = 1/14)

    container <- list(
      state = c(VH1 = 0, VH2 = 0, SH = 10000, EH = 0, IH = 5, RH = 0,
                SA = 175000, IA = 100, QA = 0, PA = 0,
                SB = 175000, IB = 100, QB = 0, PB = 0,
                VM1 = 0, VM2 = 0, SM = 100000, EM = 0, IM = 200, RM = 0),

      parameters = c(human_point_estimates,
                     livestock_point_estimates,
                     culex_point_estimates,
                     aedes_point_estimates),

      times = seq(1, 500))

      model_out <- with(container,
                          model_run <- ode(y = state, times = times,
                                           func = rvf_vaxx, parms = parameters))

      model_out <- data.frame(model_out)

      #NOTE: total humans is hardcoded
      infected_at_peak = max(model_out$IH)
      prop_infected = infected_at_peak/10005
      day_of_peak = which(model_out$IH == infected_at_peak)
      row = which(human_range == nuH)
      column = which(livestock_range == nuM)
      max_infected[row, column] <- prop_infected
      time_to_peak[row, column] <- day_of_peak
  }
}

# NOTE: my results are currently max number of infected in a single day
# At the epidemic's peak, XX% of the population is infected.

# 7/28: The value 0.00019996 is 2/10002 so no epidemic occurred
max_infected
time_to_peak # oh! that might be the more effective value. Also, no amount of
             # running this all the way through will be much more interesting

# appears to take off over a range of cattle but not over humans
# 7/28: reran with human range only over 0-10%. Still don't see an epidemic
#     past 0% human vaxx. Also see a weirdly delayed epidemic at 40% livestock
#     vaxx
