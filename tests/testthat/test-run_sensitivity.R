test_that("generate LHS sampling", {
  init = c(SH = 1000, IH = 0, RH = 0, VH = 100,
           SM = 10000, IM = 20, RM = 0, VM = 100,
           SA = 100000, IA = 200, PA = 0, QA = 0)

  parameter_names = c("nuH", "omega", "alphaH", "nuM", "muA", "alphaM", "tau",
                      "rhoRM", "betaM", "rhoRH", "rhoVM", "muE", "betaAM",
                      "betaMA", "betaHA", "rhoVH", "zeta", "betaAH")

  test <- run_sensitivity(state = init, parameter_names = parameter_names,
                          numLHS = 5000, times = 1:50, model = rvf_sensitivity,
                          rank = TRUE, nboot = 100, conf = 0.95)
  expect_equal(ncol(LHS), length(parameters))
  expect_equal(nrow(LHS), numLHS)
})

test_that("run LHS through disease model", {
  init = c(SH = 1000, IH = 0, RH = 0, VH = 100,
           SM = 10000, IM = 20, RM = 0, VM = 100,
           SA = 100000, IA = 200, PA = 0, QA = 0)

  parameter_names = c("nuH", "omega", "alphaH", "nuM", "muA", "alphaM", "tau",
                      "rhoRM", "betaM", "rhoRH", "rhoVM", "muE", "betaAM",
                      "betaMA", "betaHA", "rhoVH", "zeta", "betaAH")

  test <- run_sensitivity(state = init, parameter_names = parameter_names,
                          numLHS = 5000, times = 1:50, model = rvf_sensitivity,
                          rank = TRUE, nboot = 100, conf = 0.95)
  y <- test[[2]]

  expect_equal(length(y), numLHS)
})
