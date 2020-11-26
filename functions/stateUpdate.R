## STATE UPDATE FUNCTION

## A function which tracks the state change of units within the population over X generations. The function takes as input:  
## This requires functions S_to_I and R_to_S


library(here)
here()

source(here("functions", "effcChecker.R"))
source(here("functions", "effcIdGenerator.R"))
source(here("functions", "S_to_I.R"))


state_update <- function(state_tracker, # vector of initial population states
                         generations, # number of timesteps to run simulation 
                         R_period, # the immune period (number of timesteps spent in recovered state)
                         N, # population size
                         heterogeneity, # is heterogeneity in effective contacts "fixed" or "variable"
                         k_estimate, # overdispersion parameter, level of heterogeneity in effective contact potential
                         R_estimate, # average number of effective contacts per unit
                         output # defines the output as "matrix" or "counts"
) {
  
  
  # -------------
  ## SET UP ##
  # -------------
  
  # Create a vector to store the number of timesteps spent in R state
  # this will ensure that units spend the number of generations defined by R_period in the "R" state before reverting to susceptibility
  RTime <- c(rep(0, N))
  
  
  # Define the effective contact potential of units, when heterogeneity == "fixed"
  # If the number of effective contacts per unit is fixed over time then define a vector "effc_fixed" to contain the number of contacts per unit in the population.
  fixed_effc <- rnbinom(N,
                        size = k_estimate, 
                        mu = R_estimate)
  
  # Use effcChecker function to validate that the number of effective contacts per unit does not exceed the total population size
  fixed_effc <- effcChecker(fixed_effc, N, k_estimate, R_estimate)
  
  
  # Set the initial state of the population, this will be updated with each iteration of the generations FOR loop below.
  
  disease_state <- state_tracker
  
  
  # Create vectors to store the number of units in each state per generation, only if the output is "counts"
  
  if(output == "counts"){
    S_counts <- sum(disease_state == "S")
    I_counts <- sum(disease_state == "I")
    R_counts <- sum(disease_state == "R")
  }
  
  # --------------
  ## SIMULATION ##
  # --------------
  
  for (i in 2:generations) { # run the simulation for the number of timsteps specified in 'generations' variable
    
    S_index <- which(disease_state == "S") # ids of Susceptible individuals
    I_index <- which(disease_state == "I") # ids of Infected individuals
    R_index <- which(disease_state == "R") # ids of Recovered individuals
    
    new_disease_state <-
      vector(length = length(disease_state)) # vector to hold disease state of units in population at next timestep
    
    # ------------
    ## RECOVERY ##
    # ------------
    
    # All infecteds become recovered in a given timestep
    
    if (length(I_index) > 0) {
      new_disease_state[I_index] <- "R"
    }
    
    
    # -------------------------------
    ## LOSS OF IMMUNITY ##
    # -------------------------------
    
    ## Recovereds become Susceptible again after being in R_state for >= R_period timesteps
    
    if (length(R_index) > 0) {
      RTime[R_index] <- RTime[R_index] + 1 # update R units in RTime vector to reflect +1 timestep in R state
      newS_index <- which(RTime == R_period) # identify units which have been in R state for >= ~R_period, store in newS_index
      
      # Update RTime vector for units reverting to susceptibility
      RTime[newS_index] <- 0 # becomes 0 since units are no longer immune
      
      # Update new_disease_state vector
      new_disease_state[R_index] <- "R"
      new_disease_state[newS_index] <- "S" # units identified in newS_index lose immunity
    }
    
    
    # -----------
    ## INFECTION
    # -----------
    
    ## Infected units produce secondary cases according to negative binomial distribution (NBD)
    # See the S_to_I function description above for a detailed explanation
    # S_to_I takes a vector of susceptible units and a vector of infected units and generated new infected units (cases) based on anegative binomial distribution. 
    # S_to_I outputs a vector (cases_index) storing the ID of new infected units (cases)
    
    if (length(S_index > 0)) {
      cases_index <- S_to_I(I_index,
                            S_index,
                            N,
                            k_estimate,
                            R_estimate,
                            heterogeneity, # heterogeneity variable determines whether unit effective contacts are fixed or variable over time.
                            fixed_effc # stores the effective contact potential of units, used if heterogeneity == "fixed"
      ) 
      new_disease_state[S_index] <- "S" # all S units remain in S state
      new_disease_state[cases_index] <- "I" # update ids of units in cases_index to "I"
    }
    
    # -----------------------
    ## UPDATE STATE MATRIX // COUNTS ##
    # -----------------------
    
    if (output == "matrix"){
      
      state_tracker <-
        cbind(state_tracker, new_disease_state) # update column of state_matrix with new population state
      
    }else if (output == "counts"){
      
      # update vectors containing counts of units in each disease state over time
      
      S_counts <- c(S_counts, sum(new_disease_state == "S"))
      I_counts <- c(I_counts, sum(new_disease_state == "I"))
      R_counts <- c(R_counts, sum(new_disease_state == "R"))
      
    }
    
    # update the disease state of the population for the next loop
    disease_state <- new_disease_state
    
  } 
  ## END OF FOR LOOP
  
  # -----------------------
  ## DEFINE OUTPUT ##
  # -----------------------
  
  # Define the formate of the function output, as specified by the 'output' argument
  # If output is set to 'matrix' return the full matrix containing the disease state of each individual unit over time
  # If output is set to 'counts' return a dataframe tracking only the number (not the identity) of units in each disease state over time
  
  if (output == "matrix") {
    
    #rename state_tracker columns as generation number
    colnames(state_tracker) <- seq(1:generations)
    
    return(state_tracker) # return the full matrix tracking the state of each unit in the population across time
    
  } else if (output == "counts") {
    
    # bind state_count vectors into a data_frame called state_counts, with each disease state occupying a different row (S, I, R)
    state_counts <- as.data.frame(rbind(S_counts, I_counts, R_counts))
    # rename columns and rows of state_counts 
    colnames(state_counts) <- seq(1:generations)
    rownames(state_counts) <- c("S", "I", "R")
    
    return(state_counts) # return the number (but not ID) of units in each state across time
  }
  
}




