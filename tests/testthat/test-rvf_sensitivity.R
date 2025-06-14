test_that("no negative compartments", {

  state = c(SH = 1000, IH = 0, RH = 0, VH = 100,
            SM = 10000, IM = 20, RM = 0, VM = 100,
            SA = 100000, IA = 200, PA = 0, QA = 0)

  parameters = c(nuH = 0.01, omega = 0.2, alphaH = 1/7, nuM = 0.01, muA = 0.001,
                 alphaM = 1/7, tau = 0.2, rhoRM = 1/(365*4), betaM = 0.001,
                 rhoRH = 1/(365*7), rhoVM = 1/(365*2), muE = 0.05,
                 betaAM = 0.001, betaMA = 0.89, betaHA = 0.89,
                 rhoVH = 1/(365*2), zeta = 0.2, betaAH = 0.001)

  times = seq(1,10)

  test <- ode(y = state, times = times, func = rvf_sensitivity, parms = parameters)
  expect_true(all(test >= 0))

})




