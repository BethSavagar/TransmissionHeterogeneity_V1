

# Create a function which takes as input a vector containing the states (S, I, R) of each individual in the population  and produces a new vector containing the state (S, I, R) of each individual in the next generation

# See RF_model_draft.Rmd file for full explanation



## Function inputs: ##
# state_vector: A vector containing the disease state for each individual in a population
# p_transmission: The probability of transmission

change_state <- function(state_vector, 
                         p_transmission){
  
  # Calculation of r, the risk of a susceptible becoming infected in a given generation
  # r depends on the probability of transmission (p_transmission) and number of infectives in the current generation (seed)
  
  seed <- length(state_vector[state_vector == "I"]) # the number of infected individuals in the current generation
  r <- 1 - (1 - p_transmission)^seed # the probability of a susceptible individual becoming infected in a given generation
  
  new_state_vec <- vector(length = length(state_vector)) # A vector to store the new disease states for all individuals in the population
  
  for (i in 1:length(state_vector)) { # iterating over each individual in the population
    
    # IF an individual's state is S then run simulation to determine whether the individual becomes infected by a Bernouilli Trial (binomial distribution) with probability of infection, r. 
    
    if (state_vector[i] == "S") { 
      new_state <- rbinom(1, 1, r) # r is the likelihood of becoming infected, size = 1 (since 1 individual is being considered), n = 1 (since 1 trial)
      # the new_state variable codes an individual's state as 0 (if unchanged, S) and 1 (if changed, I)
      
    } else if (state_vector[i] == "I" | state_vector[i] == "R") { # If an individual's state is either I or R then the new state will be R
      new_state <- 2 # the new_state variable codes an individual's state as 2, denoting R
      
    }
    # Update the new_state_vec vector to contain the state of all individuals in the population in the new generation, coded as 0,1,2 for S,I,R states respectively.
    
    new_state_vec[i] <- new_state 
  }
  
  # Recode the new_state_vec vector such that 0,1,2 are replaced with strings S, I, R, as in the input state_vector, for future useability (see example below)
  
  new_state_vec[new_state_vec == 0] <- "S"
  new_state_vec[new_state_vec == 1] <- "I"
  new_state_vec[new_state_vec == 2] <- "R"
  
  return(new_state_vec) # output is the new_state_vec which contains the state of all individuals in the new generation
}

