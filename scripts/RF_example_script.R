# An example using the RF_model function (called change_state) found in functions folder
# In the example the disease state of individuals in a population of 100 is simulated over 20 generations, given a single initial index case.


library(here)
here() #setting location as .Rproj folder
source(here("functions/RF_model_fn_140820.R"))


N <- 100 # Population size
I0 <- 1 # Index case

init_state <- c(rep("I", length(I0)), 
                rep("S", N - length(I0))
) # A vector containing the state of all individuals in the populationL: a single index case (I0) and the remaining all susceptible

p_transmission <- 0.5 # 50% probability of transmission - high to check functionality!


# Test 1, given an initial state vector, the change state function will produce a vector containing the state of all individuals in the next generation

state2 <- change_state(init_state, p_transmission)
head(cbind(init_state, state2))


# Test 2, simulate the change in state of all individuals in the population over 20 generations, create a matrix to contain the state change.

generations <- 20 # number of generations to simulate

state_tracker <- matrix(nrow = N, ncol = 1+generations) # a matrix to contain the disease state of all individuals over 20 generations, ncol = 1 + generations with first col containing initial states

state_tracker[,1] <- init_state # first col contains initial state of individuals in population

# A for loop to fill the state_tracker matrix
for(i in 2:ncol(state_tracker)){ # iterate over each generation of individuals
  
  current_state <- state_tracker[,i-1] # update the current state vector to reflect the state in the current generation 
  I0 <- length(current_state[current_state == "I"]) # the seed is the number of infectives in the current generation (held in current_state)
  new_state <- change_state(current_state, p_transmission) # use the change_state function to update the state of the individuals in the population, new_state is a vector containing disease states of all individuals in the new generation
  state_tracker[, i] <- new_state # update the state_tracker matrix with the new_state vector for the new generation
  
}

head(state_tracker)
