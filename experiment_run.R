#' Run simple exploratory experiments on two-host/two-vector Rift Valley Fever 
#'  model run inputs
#'  
#' `experiment_run()` allows for exploratory experimentation by creating its own 
#'  "container" of default initial conditions, parameter values, time steps, and 
#'  function that can be modified to understand their effects when changed. 
#'  
#' @param ID A numeric value to tag output results, specified in double 
#'    quotations
#' @param variable_name Variable or concatenated variable pair to be modified 
#'    during experiment, specified by name within double quotations. Variables 
#'    that can be modified are `"state"`, `"parameters"`, `"times"`, or `"func"`
#' @param value_name Value or concatenated value pair within `variable_name` to
#'    be modified, specified by name within double quotations. Values that can 
#'    be modified are any initial state, parameter value, time step sequence, 
#'    or named Rift Valley Fever function
#' @param new_value Numeric variable value to replace default value   
#' @param brief_description A brief verbal description of experiment performed
#'    input as free text and returned with function return 
#' @returns Experiment information, including: experiment ID, changes as 
#'    specified by function inputs, and brief description provided by user. Also
#'    returns a four-panel plot showing disease behavior for each of two host
#'    and two vector species. 
#' @examples
#' library(deSolve)
#' experiment_run("001", "state", "SH", 7500, "Change initial susceptible human
#'   population from 10,000 to 7,500")
#'   
#' experiment_run("002", c("state", "parameters"), c("SH", "nuTH"), c(7500, 0.90), 
#'   "Change initial susceptible human population from 10,000 to 7,500 and 
#'   change human vaccine take from 100% to 90% effective rate")

experiment_run <- function(ID, variable_name, value_name, new_value, brief_description) {
  
  # embed default values into a "container"
  container <- list(
    state = c(VH = 0, SH = 10000, EH = 0, IH = 2, RH = 0,
              SA = 175000, IA = 100, QA = 0, PA = 0,
              SB = 175000, IB = 100, QB = 0, PB = 0,
              VM = 0, SM = 100000, EM = 0, IM = 200, RM = 0),
    
    parameters = c(nuH = 0.001, nuTH = 1, betaHA = 0.89, betaHB = 0.89, xiH = 1/4, alphaH = 1/(3*0.99),
                   rhoRH = 1/900, rhoVH = 1/360, muH = 4/(2*50*360), gammaH = 4/(2*50*360),
                   
                   nuM = 0.001, nuTM = 1, betaMA = 0.89, betaMB = 0.89, xiM = 24/3.25, alphaM = 1/(3*0.95),
                   rhoRM = 1/900, rhoVM = 1/360, muM = 0.0008,
                   gamma1 = 0.00082, gamma2 = 0.000082,
                   
                   betaAH = 0.01, betaAM = 0.01, tauA = 0.2, zetaA = 0.5, omegaA = 0.01,
                   kappaA = 175000, gammaA = 10, omegaA11 = 0.01,
                   muPA = 0.00001, muSA = 1/3, muIA = 1/3, muQA = 0.00001,
                   
                   betaBH = 0.01, betaBM = 0.01, tauB = 0.2, zetaB = 0.05, omegaB = 0.01,
                   kappaB = 175000, gammaB = 25, omegaB11 = 0.01,
                   muPB = 0.005, muSB = 0.1, muIB = 0.1, muQB = 0.005),
    
    times = seq(1, 500),
    
    func = rvf_vaxx)
  
  # edit container values based on inputs from function call:
  if(length(variable_name) == 1) {
    container2 <- container
    vector <- container2[[variable_name]]
    vector |> data.frame()
    vector[[value_name]] <- new_value
    container2[[variable_name]] <- vector
  }
  
  if(length(variable_name) == 2) {
    container1 <- container
    idx <- which(names(container1) == variable_name)
    vector <- container1[[idx[1]]]
    vector |> data.frame()
    vector[[value_name[1]]] <- new_value[1]
    container1[[variable_name[1]]] <- vector
    
    container2 <- container1
    vector <- container2[[idx[2]]]
    vector |> data.frame()
    vector[[value_name[2]]] <- new_value[2]
    container2[[variable_name[2]]] <- vector
    container2[[variable_name[2]]] <- vector
  }
  
  # run model:
  model_out <- with(container2,
                    model_run <- ode(y = state, times = times, func = rvf_vaxx, parms = parameters)
  )
  
  # bespoke plot function:
  #-- consists of two pieces: plot_host, plot_vector
  #-- each are called twice in experiment_plot()
  plot_host <- function(model_out, host) {
    
    model <- data.frame(model_out)
    
    if(host == "Human") {
      S = model$SH
      E = model$EH
      I = model$IH
      R = model$VH
      V = model$VH
      
      par(mar = c(2,4,2,4))
      plot <- plot(S ~ times, type = "l", col = "green", lwd = 2, ylab = NA,
                   xlab = "Days", ylim = c(0, S[1]+100))
      lines(E ~ times, col = "orange", lwd = 2)
      lines(I ~ times, col = "red", lwd = 2)
      lines(R ~ times, col = "purple", lwd = 2)
      lines(V ~ times, col = "turquoise", lwd = 2)
      
      legend(x = 3/4*max(times), y = S[1] - 5,
             legend=c("S","E", "I", "R", "V"),
             lwd = 4,
             col = c("green", "orange", "red", "purple", "turquoise"),
             title= host, bty = "n")
      
    }
    
    if(host == "Mammal") {
      S = model$SM
      E = model$EM
      I = model$IM
      R = model$RM
      V = model$VM
      
      par(mar = c(2,4,2,4))
      plot <- plot(S ~ times, type = "l", col = "green", lwd = 2,
                   ylab = NA, xlab = "Days", ylim = c(0, S[1]+100))
      lines(E ~ times, col = "orange", lwd = 2)
      lines(I ~ times, col = "red", lwd = 2)
      lines(R ~ times, col = "purple", lwd = 2)
      lines(V ~ times, col = "turquoise", lwd = 2)
      
      legend(x = 3/4*max(times), y = S[1] - 5,
             legend=c("S","E", "I", "R", "V"),
             lwd = 4,
             col = c("green", "orange", "red", "purple", "turquoise"),
             title= host, bty = "n")
      
    }
  }
  
  plot_vector <- function(model_out, vector) {
    
    model <- data.frame(model_out)
    
    if(vector == "A") {
      P = model$PA
      S = model$SA
      I = model$IA
      Q = model$QA
      
      par(mar = c(2,4,2,4))
      plot <- plot(S ~ times, type = "l", col = "blue", lwd = 2, xlab = "Days",
                   ylab = NA, ylim = c(0, max(c(S[1],I[1]))) + 1000)
      lines(P ~ times, col = "lightblue", lwd = 2)
      lines(I ~ times, col = "red", lwd = 2)
      lines(Q ~ times, col = "pink", lwd = 2)
      
      legend(x = 3/4*max(times), y = max(c(S[1],I[1])) - 100,
             legend=c("P","S", "I", "Q"),
             lwd = 4,
             col = c("blue", "lightblue", "red", "pink"),
             title= paste("Vector ", vector, sep = ""), bty = "n")
      
    }
    
    if(vector == "B") {
      P = model$PB
      S = model$SB
      I = model$IB
      Q = model$QB
      
      
      par(mar = c(2,4,2,4))
      plot <- plot(S ~ times, type = "l", col = "blue", lwd = 2, xlab = "Days",
                   ylab = NA, ylim = c(0, max(c(S[1],I[1]))) + 1000)
      lines(P ~ times, col = "lightblue", lwd = 2)
      lines(I ~ times, col = "red", lwd = 2)
      lines(Q ~ times, col = "pink", lwd = 2)
      
      legend(x = 3/4*max(times), y = max(c(S[1],I[1])) - 100,
             legend=c("P","S", "I", "Q"),
             lwd = 4,
             col = c("blue", "lightblue", "red", "pink"),
             title= paste("Vector ", vector, sep = ""), bty = "n")
      
    }
  }
  
  experiment_plot <- function(model) {
    
    plot_host(model, "Human")
    plot_vector(model, "A")
    plot_vector(model, "B")
    plot_host(model, "Mammal")
    
  }
  
  # autoplot output:
  par(mfrow = c(4,1))
  experiment_plot(model_out)
  
  # also return model metadata (I'd love if this actually returned first)
  return(c(paste("Experiment ID = ", ID, sep = ""),
           paste("Change ", variable_name, "$", value_name, " to value = ", new_value, sep = ""),
           brief_description))
}
