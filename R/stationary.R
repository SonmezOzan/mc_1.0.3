# stationary() finds the stationary distribution for the given markov
# chain (defined through a probability transition matrix or function)

# Input variables:
#    pijdef:  The transition probabilities, either in matrix form or a function
#    type:  Type of markov chain, either 'discrete' or 'continuous'
#    tol:  A positive scalar for error tolerance for infinite markov chain approximation
#    qidef: Holding rates at each state for continuous markov chain
#    transrate: Instead of inputting the individual holding rate and transition probabilities,
#               user can input the transition rate instead for continuous time markov chain,
#               as a function or a matrix

stationary = function(pijdef=NULL, type, tol = 1e-6, ...){
  args = list(...)
  if(type == 'discrete'){ #Discrete case
    states.type = class(pijdef)
    if(states.type == 'matrix'){ #Finite number of states
      mkc = new('markovchain', transitionMatrix = pijdef)
      absorb = absorbingStates(mkc)
      if(length(absorb) > 0) stop('At least 1 absorbing state(s). Use absorb.mc() instead.')
      
      pis = findpil.fin(pijdef)
      mc = list(pijdef = pijdef, stationary.distribution = pis)
      class(mc) = 'mc'
      return(mc)
    }else{ #Infinte number of states
      pis = findpil.inf(pijdef, tol)
      mc = list(pijdef = pijdef, stationary.distribution = pis)
      class(mc) = 'mc'
      return(mc)
    }
    
  }else if(type == 'continuous'){ #Continuous case. 
    if(!hasArg('qidef') && !hasArg(transrate)) stop('Missing holding/transition rates')
    if(hasArg(transrate)){ #User input transition rate
      state.type = class(args$transrate)
      if(state.type == 'matrix'){ #Finite states
        pis = findpicont.fin(transrate = args$transrate)
        mc = list(transrate = args$transrate, stationary.distribution = pis)
        class(mc) = 'mc'
        return(mc)
      }else{ #Infinite states
        pis = findpicont.inf(transrate = args$transrate)
        mc = list(transrate = args$transrate, stationary.distribution = pis)
        class(mc) = 'mc'
        return(mc)
      }
    }else{ #User input probability matrix/function and holding rates
      if(is.null(pijdef)) stop('Missing probability matrix')
      state.type = class(pijdef)
      if(state.type == 'matrix'){ #Finite number of states
        if(length(args$qidef) != nrow(pijdef)) stop('Dimension of probability matrix and holding rates mismatch')
        pis = findpicont.fin(pijdef = pijdef, qidef = args$qidef)
        mc = list(pijdef = pijdef, stationary.distribution = pis, holding.rates = args$qidef)
        class(mc) = 'mc'
        return(mc)
      }else{
        pis = findpicont.inf(pijdef = pijdef, qidef = args$qidef)
        mc = list(pijdef = pijdef, stationary.distribution = pis, holding.rates = args$qidef)
        class(mc) = 'mc'
        return(mc)
      }
    }
  }
}


#Find the stationary vector given a finite transition probability matrix
findpil.fin = function(pijdef){ 
  n <- nrow(pijdef)
  imp <- diag(n) - t(pijdef)
  imp[n, ] <- rep(1, n)
  rhs <- c(rep(0, n-1), 1)
  solve(imp, rhs)
}

#Find the stationary probability given a function for infinite state transition probability
findpil.inf = function(pijdef, tol = 1e-06){
  k = 10
  pij = sapply(1:k, function(i){
    sapply(1:k, function(j){
      pijdef(i,j)
    })
  })
  pij = t(pij)
  pij[k,k] = 1-sum(pij[k,1:(k-1)])
  stationary.pi = findpil.fin(pij)
  if(stationary.pi[length(stationary.pi)] > stationary.pi[length(stationary.pi)-1]) stop("stationary distribution doesn't converge")
  if(all(abs(diff(stationary.pi)) < 1e-10)) stop('no stationary distribution')
  
  k = 20
  pij = sapply(1:k, function(i){
    sapply(1:k, function(j){
      pijdef(i,j)
    })
  })
  pij = t(pij)
  pij[k,k] = 1-sum(pij[k,1:(k-1)])
  stationary.pi.new = findpil.fin(pij)
  error = sum((c(stationary.pi, rep(0, 10)) - stationary.pi.new)^2)
  
  while(error > tol){
    stationary.pi = stationary.pi.new
    k = k + 10
    pij = sapply(1:k, function(i){
      sapply(1:k, function(j){
        pijdef(i,j)
      })
    })
    pij = t(pij)
    pij[k,k] = 1-sum(pij[k,1:(k-1)])
    stationary.pi.new = findpil.fin(pij)
    error = sum((c(stationary.pi, rep(0, 10)) - stationary.pi.new)^2)
  }
  return(stationary.pi.new)
}

#Find the stationary distribution for continuous markov chain, given the transition probability and holding rates
findpicont.fin = function(...){
  args = list(...)
  if(hasArg(pijdef)){ #User input probability matrix and holding rate
    #First create the Q matrix
    n = length(args$qidef)
    Q = diag(-args$qidef)
    for(i in 1:n){
      Q[-i, i] = args$qidef[i]*args$pijdef[i,-i]
    }
    Q[n, ] = rep(1,n)
    rhs = c(rep(0, n-1), 1)
    pivec = solve(Q, rhs)
    return(pivec)
  }else{
    Q = args$transrate
    n = nrow(Q)
    Q[n,] = rep(1,n)
    rhs = as.matrix((c(rep(0,n-1),1)))
    pivec = solve(Q,rhs)
    return(pivec)
  }
}

#Find the stationary distribution for infinite state continuous markov chain, given the transition probability and holding rates
findpicont.inf = function(tol = 1e-6,...){
  args = list(...)
  
  truncate = function(k, transrate){
    Q = sapply(1:k, function(i){
      sapply(1:k, function(j){
        transrate(i,j)
      })
    })
    Q = t(Q)
    return(Q)
  }
  
  if(hasArg(transrate)){ #User input transition rates
    Q = truncate(10, args$transrate)
    stationary.old = findpicont.fin(transrate = Q)
    Q = truncate(20, args$transrate)
    stationary.new = findpicont.fin(transrate = Q)
    error = sum((c(stationary.old, rep(0, 10)) - stationary.new)^2)
    k = 20
    if(stationary.new[length(stationary.new)] > stationary.new[length(stationary.new)-1]) stop("stationary distribution doesn't converge")
    if(all(abs(diff(stationary.pi)) < 1e-10)) stop('no stationary distribution')
    
    while(error > tol){
      stationary.old = stationary.new
      k = k + 10
      Q = truncate(k, args$transrate)
      stationary.new = findpicont.fin(transrate = Q)
      error = sum((c(stationary.old, rep(0, 10)) - stationary.new)^2)
    }
    return(stationary.new)
  }else{ #User input functions for probability matrix and holding rates
    rate.truncate = function(k, rates){
      qi = sapply(1:k, function(i){
        rates(i)
      })
      return(qi)
    }
    
    pij = truncate(10, args$pijdef)
    qi = rate.truncate(10, args$qidef)
    stationary.old = findpicont.fin(pijdef = pij, qidef = qi)
    pij = truncate(20, args$pijdef)
    qi = rate.truncate(20, args$qidef)
    stationary.new = findpicont.fin(pijdef = pij, qidef = qi)
    error = sum((c(stationary.old, rep(0, 10)) - stationary.new)^2)
    k = 20
    
    while(error > tol){
      stationary.old = stationary.new
      k = k + 10
      pij = truncate(k, args$pijdef)
      qi = rate.truncate(k, args$qidef)
      stationary.new = findpicont.fin(pijdef = pij, qidef = qi)
      error = sum((c(stationary.old, rep(0, 10)) - stationary.new)^2)
    }
    return(stationary.new)
  }
}