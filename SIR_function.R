
# SIR Epidemic model function:

SIR_dyn <- function(time, states, params){
  
  # t <- 
  # states <- initial state vector for S, I, R compartments
  # par <- vector of parameters, beta and gamma
  
  with(as.list(c(states, params)), {
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  })
  
}
