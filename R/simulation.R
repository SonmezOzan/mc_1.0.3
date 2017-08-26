
# Discrete Markov Chain Simulation for a given transition matrix 
# and the number of desired stages of the markov chain

# Input variables
   # pijdf: The transition probabilities, either in matrix form or a function
   #type:  Type of markov chain, either 'discrete' or 'continuous'
   # N = Total number of steps that MC will run
   # initial.state : initial state that the chain will start from, if not specified then
   #it will be chosen randomly from the state space.

mc.simulation <- function(pijdef, type, N, initial.state=NULL, ...){
   args = list(...)
   if(type == 'discrete'){ #Discrete case
      states.type = class(pijdef)
      if(states.type == 'matrix'){ #Finite number of states
         mkc = new('markovchain', transitionMatrix = pijdef)
         colnames(pijdef) <- c(1:ncol(pijdef))
         rownames(pijdef) <- c(1:nrow(pijdef))
         
         transit <- function(char, pijdef) {
            sample(colnames(pijdef), 1, prob = pijdef[char,])
         }
         sim <- character(N)
         if (is.null(initial.state)){
            sim[1] <- sample(colnames(pijdef), 1)
         }else{
            sim[1] <- initial.state
         }
         for (i in 2:N) {
            sim[i] <- transit(sim[i-1], pijdef)
         }
         simulation = list(pijdef = pijdef, simulated.states = as.numeric(sim))
         return(simulation)
      }else{
         print('infinite state space (to be handled)')
      }
      
   }else if(type == 'continuous'){ #Continuous case. 
      print('to be continued')
   }
}
