build_state <- function(SH = 10000, IH = 2, 
                        SM = 100000, IM = 200) {
  #' @param SH Susceptible humans as count 
  #' @param IH Infected humans as count 
  #' @param SM Susceptible livestock as count
  #' @param IM Infected livestock as count
  #' @returns String containing initial state values as proportions
  #' for input into freqency-dependent SIR model
  
  # Input initial states as intuitive counts
  human_count = c(VH1 = 0, VH2 = 0, SH = SH, EH = 0, IH = IH, RH = 0)
  livestock_count = c(VM1 = 0, VM2 = 0, SM = SM, EM = 0, IM = IM, RM = 0)
  culex_count = c(SA = 175000, IA = 100, QA = 0, PA = 0)
  aedes_count = c(SB = 175000, IB = 100, QB = 0, PB = 0)
  
  # Convert to densities = count/N
  human_density = human_count/sum(human_count)
  livestock_density = livestock_count/sum(livestock_count)
  culex_density = culex_count/sum(culex_count)
  aedes_density = aedes_count/sum(aedes_count)
  
  # Must be in model output order
  state = c(human_density,
            culex_density,
            aedes_density,
            livestock_density)
  
  return(state)
}
 
build_parameters <- function(nuH = 0, nuM = 0) {
  
  #' @param nuH Human vaccination rate in percentage 
  #' @param nuM Livestock vaccination rate in percentage 
  #' @returns String containing default parameter values
  
  parameters = c(nuH = nuH, nuTH = 1, betaHA = 0.89, betaHB = 0.89,
                 xiH = 1/4, alphaH = 1/(3*0.99), rhoRH = 1/900, rhoVH = 1/360,
                 muH = 4/(2*50*360), gammaH = 4/(2*50*360), thetaH = 1/14,
  
                 nuM = nuM, nuTM = 1, betaMA = 0.89, betaMB = 0.89,
                 xiM = 24/3.25, alphaM = 1/(3*0.95), rhoRM = 1/900,
                 rhoVM = 1/360, muM = 0.0008, gamma1 = 0.00082,
                 gamma2 = 0.000082, thetaM = 1/14,
                 
                 betaAH = 0.01, betaAM = 0.01, tauA = 0.2, zetaA = 0.5,
                 omegaA = 0.01, kappaA = 175000, gammaA = 10, omegaA11 = 0.01,
                 muPA = 0.00001, muSA = 1/3, muIA = 1/3, muQA = 0.00001,
                 
                 betaBH = 0.01, betaBM = 0.01, tauB = 0.2, zetaB = 0.05,
                 omegaB = 0.01, kappaB = 175000, gammaB = 25, omegaB11 = 0.01,
                 muPB = 0.005, muSB = 0.1, muIB = 0.1, muQB = 0.005)
  
  return(parameters)
} 

state = build_state()
parameters = build_parameters()

