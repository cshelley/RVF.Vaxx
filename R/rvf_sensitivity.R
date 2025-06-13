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
#' parameters = c(nuH = 0.1*0.9, rhoVH = 1/(365*7), rhoRH = 1/(365*10),
#'                betaAH = 0.001, alphaH = 1/7,
#'                nuM = 0.1*0.9, rhoVM = 1/(365*2), rhoRM = 1/(365*10),
#'                betaAM = 0.001, alphaM = 1/7,
#'                tau = 0.2, omega = 0.01, zeta = 0.05, betaHA = 0.89,
#'                betaMA = 0.89, muA = 0.1, muE = 0.005, kappaA = 175000)
#'
#' times = seq(1, 500)
#'
#' # run model and convert output to a data.frame
#' ode(y = state, times = times, func = rvf_sensitivity, parms = parameters) |>
#'   data.frame()

rvf_sensitivity <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{

    # support equations
    NH = VH + SH + IH + RH
    NM = VM + SM + IM + RM
    Na = SA + IA

    # Human Hosts:
    dVH = nuH*(SH + IH + RH) - rhoVH*VH
    dSH = rhoVH*VH - betaAH*SH*IA/Na - nuH*SH + rhoVH*VH + rhoRH*RH
    dIH = betaAH*SH*IA/Na - alphaH*IH - nuH*IH
    dRH = alphaH*IH - nuH*RH - rhoRH*RH

    # Mammalian Hosts:
    dVM = nuM*(SM + IM + RM) - rhoVM*VM
    dSM = rhoVM*VM - betaAM*SM*IA/Na - nuM*SM + rhoVM*VM + rhoRM*RM
    dIM = betaAM*SM*IA/Na - alphaM*IM - nuM*IM
    dRM = alphaM*IM - nuM*RM - rhoRM*RM

    # Vector:
    dSA = -betaHA*SA*IH/NH - betaMA*SA*IM/NM + tau*PA - muA*SA
    dIA = betaHA*SA*IH/NH + betaMA*SA*IM/NM + tau*QA - muA*IA
    dPA = omega*SA + (1-zeta)*omega*(1-Na/kappaA)*IA - tau*PA - muE*PA
    dQA = zeta*omega*(1-Na/kappaA)*IA - tau*QA - muE*QA

    # return rates of change
    list(c(dVH, dSH, dIH, dRH,
           dSA, dIA, dQA, dPA,
           dVM, dSM, dIM, dRM))
  })
}
