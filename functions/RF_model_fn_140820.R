


# Function 1: given state vector transform S --> I

change_state <- function(state_matrix,
                         state_vector, 
                         p_transmission){
  
  new_state_vec <- vector(length = length(state_vector)) # A vector to store the new disease states for all individuals in the population
  
  
  # Calculation of r, the risk of a susceptible becoming infected in a given generation
  # r depends on the probability of transmission (p_transmission) and number of infectives in the current generation (seed)
  
  seed <- sum(state_vector == "I") # the number of infected individuals in the current generation
  r <- 1 - (1 - p_transmission)^seed # the probability of a susceptible individual becoming infected in a given generation
  
  # IF an individual's state is S then run simulation to determine whether the individual becomes infected by a Bernouilli Trial (binomial distribution) with probability of infection, r. 
  
  S <- which(state_vector == "S") # Selects index of all individuals with state S
  
  if(length(S)>0){ # If number of S individuals is non-zero
    new_state_vec[S] <- rbinom(length(S), 1, r) # r is the likelihood of becoming infected, number of observations (n) = number of S state individuals, number of trials (size) = 1
  }
  
  # If an individual's state is I then replace with R (coded 2)
  
  new_state_vec[state_vector == "I" | state_vector == "R" ] <- 2
  
  # If an individual's state is R and has been R for >3 generations then replace with S (0)
  
  # Write a function which counts the number of Rs in a vector
  R_counter <- function(x){
    length(which(x == "R"))
  }
  
  # Use the apply family to apply the R_counter function over all rows of the state_matrix
  R_counts <- apply(state_matrix, MARGIN = 1, R_counter)
  R2S <- which(R_counts > 3) # vector of indices for individuals which have been recovered for more than 3 generations (number 3 is arbitrary)
  new_state_vec[R2S] <- 0
  
  # the new_state variable codes an individual's state as 0 (if unchanged, S) and 1 (if changed, I)
  
  # Recode the new_state_vec vector such that 0,1,2 are replaced with strings S, I, R, as in the input state_vector, for future useability (see example below)
  
  new_state_vec[new_state_vec == 0] <- "S"
  new_state_vec[new_state_vec == 1] <- "I"
  new_state_vec[new_state_vec == 2] <- "R"
  
  return(new_state_vec) # output is the new_state_vec which contains the state of all individuals in the new generation
}



#If an individual's state is R then revert to S with probability... ?

#R <- which(state_vector == "R") # vector containing indices of all R individuals



