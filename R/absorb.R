# absorbingmc() function to deal with absorbing markov chains
# put matrix into canonical form
# find the fundamental matrix N = (I-Q)^-1, N is a square matrix 
#with rows and columns corresponding to the non-absorbing states
#Thus, row i corresponds to the ith non-absorbing state, not the ith state overall.
#N(i,j) is the expected number of periods that the chain spends in the jth non-absorbing 
# the sum of each row i reveals the expected number of periods spent in any non-absorbing state 
#given that the chain initially occupied the ith non-absorbing state.
#sum(N(i,j)) is the the expected number of periods before absorption (into any absorbing state) 
#given that the chain began in the ith non-absorbing state.
#state given that the chain began in the ith non-absorbing state.

#Let bij be the probability that an absorbing chain will be absorbed in the 
#absorbing state sj if it starts in the transient state si. 
#Let B be the matrix with entries bij. Then B is an t-by-r matrix, and
#B = NR

# arguments:
#    pijdef:  The transition probabilities, either in matrix form or a function
#    type:  Type of markov chain, either 'discrete' or 'continuous'
#    tol:  A positive scalar for error tolerance for infinite markov chain approximation

# return value:
#    Object of class "mc", with components pijdef, and state vectors


absorb.mc = function(pijdef, type = 'discrete', tol = 1e-6, ...){
  args = list(...)
  if(type == 'discrete'){ #Discrete case
    states.type = class(pijdef)
    if(states.type == 'matrix'){ #Finite number of states
      mkc = new('markovchain', transitionMatrix = pijdef)
      n = nrow(mkc)
      astates = absorbingStates(mkc)
      tstates = transientStates(mkc)
      abs = which(apply(pijdef,1,function(row) 1%in%row))
      trans = which(!apply(pijdef,1,function(row) 1%in%row))
      absorb = TRUE
      if(is.irreducible(mkc)){
        absorb = FALSE
        stop('Chain is irreducible and not an absorbing Markov chain')}
      if(length(unique(astates))+length(unique(tstates))<n){
        absorb = FALSE
        stop('Chain is not an absorbing Markov chain, check that there is atleast one absorbing state and that all non-absorbing states are transient')}
      canForm = pijdef[c(trans,abs),c(trans,abs)]
      Qmat = pijdef[trans,trans]
      Rmat = pijdef[trans,abs]
      Nmat = solve(diag(nrow(Qmat)) - Qmat)
      c = matrix(rep(1),nrow(Nmat))
      times = Nmat%*%c
      colnames(times) =  'MeanTimeToAbsorb'
      rownames(times) = sapply(trans,function(x){paste('State',x)})
      abprob = Nmat%*%Rmat
      colnames(abprob) = sapply(abs,function(x){paste('State',x)})
      rownames(abprob) = sapply(trans,function(x){paste('State',x)})
      mc = list(pijdef = pijdef, absorb = absorb,astates = astates, tstates = tstates,canonicalForm = canForm,steady.state = steadyStates(mkc), 
                        FundamentalMatrix = Nmat, Qmat = Qmat, Rmat = Rmat, ExpectedHits = times,
                        AbsorbProb = abprob)
      class(mc) = 'mc'
      return(mc)
    }else{ #Infinte number of states
      print('Infinite number of states, to be continued')
    }
    
  }else if(type == 'continuous'){ #Continuous case. 
    print('to be continued')
  }
}
