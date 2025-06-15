#' Simplified Rift Valley Fever model with two hosts, two vectors, and
#'  vaccination
#'
#' `rvf_sensitivity()` is a function for in the form of a [deSolve::ode()]
#' model prepared for sensitivity analysis. It is a simplified two-host/one-
#' vector model simulating Rift Valley Fever transmission. Host species progress
#' through Susceptible (S), Infected (I), and Recovered (R), along with a
#' Vaccinated (V) state. Vector species progress from Susceptible (S) to
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
#' state = c(VH = 0, SH = 10000, IH = 2, RH = 0,
#'           SA = 175000, IA = 100, QA = 0, PA = 0,
#'           VM = 0, SM = 100000, IM = 200, RM = 0)
#'
#' parameters = c(nuH = 0.01, omega = 0.2, alphaH = 1/7, nuM = 0.01, muA = 0.001,
#'                alphaM = 1/7, tau = 0.2, rhoRM = 1/(365*4), betaM = 0.001,
#'                rhoRH = 1/(365*7), rhoVM = 1/(365*2), muE = 0.05,
#'                betaAM = 0.001, betaMA = 0.89, betaHA = 0.89,
#'                rhoVH = 1/(365*2), zeta = 0.2, betaAH = 0.001)
#'
#' times = seq(1, 500)
#'
#' # run model and convert output to a data.frame
#' ode(y = state, times = times, func = rvf_sensitivity, parms = parameters) |>
#'   data.frame()

rvf_sensitivity <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {

    # Support equations:
    Na = SA + IA
    NH = SH + IH + RH + VH
    NM = SM + IM + RM + VM

    # Human host:
    dSH <- -betaAH*SH*IA/Na + rhoRH*RH + rhoVH*VH -nuH*SH
    dIH <- betaAH*SH*IA/Na - alphaH*IH - nuH*IH
    dRH <- alphaH*IH - rhoRH*RH - nuH*RH
    dVH <- -rhoVH*VH + nuH*(SH+IH+RH)

    # Mammal host:
    dSM <- -betaAM*SM*IA/Na + rhoRM*RM + rhoVM*VM -nuM*SM
    dIM <- betaAM*SM*IA/Na - alphaM*IM - nuM*IM
    dRM <- alphaM*IM - rhoRM*RM - nuM*RM
    dVM <- -rhoVM*VM + nuM*(SM+IM+RM)

    # Vector A:
    dSA <- -betaHA*SA*IH/NH -betaMA*SA*IM/NM + tau*PA - muA*SA
    dIA <- betaHA*SA*IH/NH + betaMA*SA*IM/NM + tau*QA - muA*IA
    dPA <- -tau*PA + (1-zeta)*omega*IA + omega*SA - muE*PA
    dQA <- -tau*QA + zeta*omega*IA - muE*QA

    # return rates of change
    list(c(dSH, dIH, dRH, dVH,
           dSM, dIM, dRM, dVM,
           dSA, dIA, dPA, dQA))
  })
}
