library(FME)

run_sensitivity <- function(state, parameter_names, numLHS, times, model,
                            rank = TRUE, nboot, conf) {

  # Generate Latin Hypercube Sampling
  parRange <- data.frame(min = rep(0, length(parameter_names)),
                         max = rep(1, length(parameter_names)),
                         row.names = parameter_names)

  LHS <- data.frame(Latinhyper(parRange, numLHS))

  # Run LHS through disease model
  y <- vector()
  for(i in 1:nrow(LHS)) {
    rvf <- ode(y = init, times = times, func = rvf_sensitivity, parms = LHS[i,])
    y <- append(y, rvf[50,3], after = length(y))
  }

  # Partial rank correlation calculation
  PCC <- pcc(X = LHS, y = unlist(y), rank = rank, nboot = nboot, conf = conf)

  return(list(LHS, y, PCC$PRCC))
}


