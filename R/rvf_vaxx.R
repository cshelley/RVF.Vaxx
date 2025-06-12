#' Baseline Rift Valley Fever model with two hosts, two vectors, and vaccination
#'
#' `rvf_vaxx()` is a function for input into a [deSolve::ode()] call to run a
#' deterministic compartmental model simulating Rift Valley Fever transmission
#' in a two host/two vector system of humans (H) and mammalian hoofstock (M)
#' with Culex (A) and Aedes (B) vectors. Host species progress through
#' Susceptible (S), Exposed (E), Infected (I), and Recovered (R), along with a
#' Vaccinated (V) state. Vector species progres from Susceptible (S) to
#' Infectious (I) and can lay non-infected eggs (P) or infected eggs (Q).
#'
#' @param t Time steps to run model, expressed as a sequence from 1 to last time
#'  step
#' @param state A vector of initial conditions for all species:disease state
#'  pairs
#' @param parameters A vector of parameter values express rates in and out of
#'  model states
#' @returns A list of components corresponding to each model state of length
#'  1 through max(time step)
#' @examples
#' library(deSolve)
#'
#' state = c(VH1 = 0, VH2 = 0, SH = 10000, EH = 0, IH = 2, RH = 0,
#'           SA = 175000, IA = 100, QA = 0, PA = 0,
#'           SB = 175000, IB = 100, QB = 0, PB = 0,
#'           VM1 = 0, VM2 = 0, SM = 100000, EM = 0, IM = 200, RM = 0)
#'
#' parameters = c(nuH = 0.001, nuTH = 1, betaHA = 0.89, betaHB = 0.89,
#'                xiH = 1/4, alphaH = 1/(3*0.99), rhoRH = 1/900, rhoVH = 1/360,
#'                muH = 4/(2*50*360), gammaH = 4/(2*50*360), thetaH = 1/14,
#'
#'                nuM = 0.001, nuTM = 1, betaMA = 0.89, betaMB = 0.89,
#'                xiM = 24/3.25, alphaM = 1/(3*0.95), rhoRM = 1/900,
#'                rhoVM = 1/360, muM = 0.0008, gamma1 = 0.00082,
#'                gamma2 = 0.000082, thetaM = 1/14,
#'
#'                betaAH = 0.01, betaAM = 0.01, tauA = 0.2, zetaA = 0.5,
#'                omegaA = 0.01, kappaA = 175000, gammaA = 10, omegaA11 = 0.01,
#'                muPA = 0.00001, muSA = 1/3, muIA = 1/3, muQA = 0.00001,
#'
#'                betaBH = 0.01, betaBM = 0.01, tauB = 0.2, zetaB = 0.05,
#'                omegaB = 0.01, kappaB = 175000, gammaB = 25, omegaB11 = 0.01,
#'                muPB = 0.005, muSB = 0.1, muIB = 0.1, muQB = 0.005)
#'
#' times = seq(1, 500)
#'
#' # run model and convert output to a data.frame
#' ode(y = state, times = times, func = rvf_vaxx, parms = parameters) |>
#'   data.frame()

rvf_vaxx <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{

    # support equations
    NH = VH1 + VH2 + SH + EH + IH + RH
    NM = VM1 + VM2 + SM + EM + IM + RM
    Na = SA + IA
    NB = SB + IB
    omega2 = omegaA + omegaA11

    # Human Hosts:
    dVH1 = nuH*nuTH*(SH + EH + IH + RH) - muH*VH1 - thetaH*VH1
    dVH2 = -muH*VH2 + thetaH*VH1 - rhoVH*VH2
    dSH = gammaH*NH - betaAH*SH*IA/Na - betaBH*SH*IB/NB + rhoRH*RH - nuH*nuTH*SH + thetaH*VH2 - muH*SH
    dEH = betaAH*SH*IA/Na + betaBH*SH*IB/NB - xiH*EH - nuH*nuTH*EH - muH*EH
    dIH = xiH*EH - alphaH*IH - nuH*nuTH*IH - muH*IH
    dRH = alphaH*IH - rhoRH*RH - nuH*nuTH*RH - muH*RH

    # Mammalian Hosts:
    dVM1 = nuM*nuTM*(SM + EM + IM + RM) - muM*VM1 - thetaM*VM1
    dVM2 = -muM*VM2 +thetaM*VM1 - rhoVM*VM2
    dSM = -betaAM*SM*IA/Na - betaBM*SM*IB/NB + rhoRM*RM - nuM*nuTM*SM + rhoVM*VM - muM*SM + (gamma1*(NM-IM) + gamma2*IM)
    dEM = betaAM*SM*IA/Na + betaBM*SM*IB/NB - xiM*EM - nuM*nuTM*EM - muM*EM
    dIM = xiM*EM - alphaM*IM - nuM*nuTM*IM - muM*IM
    dRM = alphaM*IM - rhoRM*RM - nuM*nuTM*RM - muM*RM

    # Vector A:
    dSA = -betaHA*SA*IH/NH - betaMA*SA*IM/NM + tauA*PA - muSA*SA
    dIA = betaHA*SA*IH/NH + betaMA*SA*IM/NM + tauA*QA - muIA*IA
    dQA = zetaA*omegaA*(1-Na/kappaA)*gammaA*IA - tauA*QA - muQA*QA
    dPA = omega2*SA - tauA*PA + (1-zetaA)*omegaA*(1-Na/kappaA)*gammaA*IA - muPA*PA

    # Vector B:
    dSB = -betaHB*SB*IH/NH - betaMB*SB*IM/NM + tauB*PB - muSB*SB
    dIB = betaHB*SB*IH/NH + betaMB*SB*IM/NM + tauB*QB - muIB*IB
    dQB = zetaB*omegaB*(1-NB/kappaB)*gammaB*IB - tauB*QB - muQB*QB
    dPB = omega2*SB - tauB*PB + (1-zetaB)*omegaB*(1-NB/kappaB)*gammaB*IB - muPB*PB


    # return rates of change
    list(c(dVH1, dVH2, dSH, dEH, dIH, dRH,
           dSA, dIA, dQA, dPA,
           dSB, dIB, dQB, dPB,
           dVM1, dVM2, dSM, dEM, dIM, dRM))
  })
}
