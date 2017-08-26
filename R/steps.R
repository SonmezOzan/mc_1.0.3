#Usage
#Input Variables:
#pijdef: either a matrix (finite states) or a function (infinite states) capturing the 
#        transition probabilities
#type: specify whether the chain is discrete or continuous
#rows: the set of initial states
#columns: the set of final states
#tol: error tolerate for inifite states markov chain
#qidef: if the markov chain is of continuous type, users can input the transition probability via
#       pijdef, and supply the holding rates at each state using qidef
#transrate: alternatively, users can supply the transition rates, with positive values for rates
#           between different states, and -qi for rates of same state

expectstep = function(pijdef, type, rows, columns, tol = 1e-6, ...){
  args = list(...)
  rows = sort(unique(rows)); columns = sort(unique(columns))
  
  if(type == 'discrete'){ #Discrete case
    states.type = class(pijdef)
    if(states.type == 'matrix'){ #Finite number of states
      #Check if the chain is irreducible
      mkc = new('markovchain', transitionMatrix = pijdef)
      if(!is.irreducible(mkc)) warning('Chain is not irreducible, expected steps is infinite for certain states')
      absorb = absorbingStates(mkc)
      if(length(absorb) > 0){
        if(any(!(columns %in% as.numeric(absorb)))) {
          stop('Expected steps is infinite for non-absorbing states. Use absorb.mc() instead.')
        }else{
          step.matrix = steps.fin(pijdef, type, rows, columns)
        }
      }else{
        step.matrix = steps.fin(pijdef, type, rows, columns)
      }
      mc = list(pijdef = pijdef, steps = step.matrix)
      class(mc) = 'mc'
      return(mc)
    }else{ #Infinte number of states
      step.matrix = steps.inf(pijdef, type, rows, columns, tol)
      mc = list(pijdef = pijdef, steps = step.matrix)
      class(mc) = 'mc'
      return(mc)
    }
    
  }else if(type == 'continuous'){ #Continuous case. 
    if(!hasArg('qidef') && !hasArg(transrate)) stop('Missing holding/transition rates')
    if(hasArg(transrate)){ #User input transition rate
      state.type = class(args$transrate)
      if(state.type == 'matrix'){ #Finite states
        pis = timefin(rows = rows, columns = columns, transrate = args$transrate)
        mc = list(transrate = args$transrate, time = pis)
        class(mc) = 'mc'
        return(mc)
      }else{ #Infinite states
        pis = timeinf(rows = rows, columns = columns, transrate = args$transrate, tol = tol)
        mc = list(transrate = args$transrate, time = pis)
        class(mc) = 'mc'
        return(mc)
      }
    }else{ #User input probability matrix/function and holding rates
      if(is.null(pijdef)) stop('Missing probability matrix')
      state.type = class(pijdef)
      if(state.type == 'matrix'){ #Finite number of states
        if(length(args$qidef) != nrow(pijdef)) stop('Dimension of probability matrix and holding rates mismatch')
        pis = timefin(rows = rows, columns = columns, pijdef = pijdef, qidef = args$qidef)
        mc = list(pijdef = pijdef, time = pis, holding.rates = args$qidef)
        class(mc) = 'mc'
        return(mc)
      }else{
        pis = timeinf(rows = rows, columns = columns, pijdef = pijdef, qidef = args$qidef, tol = tol)
        mc = list(pijdef = pijdef, time = pis, holding.rates = args$qidef)
        class(mc) = 'mc'
        return(mc)
      }
    }
  }
}

#Find the expected number of steps for finite states
steps.fin = function(pijdef, type, rows, columns){
  n = nrow(pijdef)
  if(max(rows) > n || max(columns) > n) stop('Input states must be within range.')
  if(min(rows) < 1 || min(columns) < 1) stop('Start states at 1')
  mkc = new('markovchain', transitionMatrix = pijdef)
  absorb = absorbingStates(mkc)
  if(length(absorb) > 0){
    absorb.state = as.numeric(absorb)
    steps = sapply(columns, function(j){
      tmp = pijdef[,-j]; tmp = tmp[-j,]
      steps.j = solve(diag(n-1) - tmp, rep(1, n-1))
      return(append(steps.j, 1, after = j-1))
    })
    steps = steps[rows, ]
    steps = matrix(steps, nrow = length(rows))
    rownames(steps) = rows; colnames(steps) = columns
    return(steps)
  }else{
    pis = stationary(pijdef, type)$stationary.distribution
    steps = sapply(columns, function(j){
      tmp = pijdef[,-j]; tmp = tmp[-j,]
      steps.j = solve(diag(n-1) - tmp, rep(1, n-1))
      return(append(steps.j, 1/pis[j], after = j-1))
    })
    steps = steps[rows, ]
    steps = matrix(steps, nrow = length(rows))
    rownames(steps) = rows; colnames(steps) = columns
    return(steps)
  }
}

#Find the expected number of steps for inifite states
steps.inf = function(pijdef, type, rows, columns, tol=1e-6){
  if(min(rows) < 1 || min(columns) < 1) stop('Start states at 1')
  
  #Truncate first case
  k = max(rows, columns)
  pij = sapply(1:k, function(i){
    sapply(1:k, function(j){
      pijdef(i,j)
    })
  })
  pij = t(pij)
  rs = which(rowSums(pij) != 1)
  for(i in rs){
    pij[i,k] = 1-sum(pij[i,1:(k-1)])
  }
  
  #Check if the chain will stop or not
  check = stationary(pij, type)$stationary.distribution
  if(check[length(check)] > check[length(check)-1] && any(abs(diff(check)) >= 1e-10)) {
    if(min(rows) >= min(columns)) stop('no finite steps if i <= j')
    else warning('no finite steps if i <= j')
  }
  if(all(abs(diff(check)) < 1e-10)){
    if(any(rows %in% columns)) stop('no finite steps for i = i')
    else if(any(sapply(rows, function(i) i > columns))) stop('no finite steps for i < j')
    else warning('no finite steps for i = i')
  } 
  
  steps.old = steps.fin(pij, type, rows, columns) #Pass check, move on
  
  k = k + 50
  pij = sapply(1:k, function(i){
    sapply(1:k, function(j){
      pijdef(i,j)
    })
  })
  pij = t(pij)
  rs = which(rowSums(pij) != 1)
  for(i in rs){
    pij[i,k] = 1-sum(pij[i,1:(k-1)])
  }
  #pij[k,k] = 1-sum(pij[k,1:(k-1)])
  steps.new = steps.fin(pij, type, rows, columns)
  
  error = norm(steps.old - steps.new, type = 'O')
  while(error > tol){
    steps.old = steps.new
    k = k + 50
    pij = sapply(1:k, function(i){
      sapply(1:k, function(j){
        pijdef(i,j)
      })
    })
    pij = t(pij)
    rs = which(rowSums(pij) != 1)
    for(i in rs){
      pij[i,k] = 1-sum(pij[i,1:(k-1)])
    }
    #pij[k,k] = 1-sum(pij[k,1:(k-1)])
    steps.new = steps.fin(pij, type, rows, columns)
    error = norm(steps.old - steps.new, type = 'O')
  }
  return(steps.new)
}

#Find the expected time for finite states
timefin = function(rows, columns, ...){
  args = list(...)
  if(hasArg(pijdef)){ #User input probability matrix and holding rate
    #First create the Q matrix
    n = length(args$qidef)
    Q = diag(args$qidef)
    Q = Q %*% args$pijdef
    diag(Q) = -args$qidef
    pis = stationary(pijdef = args$pijdef, type = 'continuous', qidef = args$qidef)$stationary.distribution
    times = sapply(columns, function(j){
      tmp = Q[,-j]; tmp = tmp[-j,]
      times = solve(-tmp, rep(1,n-1))
      return(append(times, 1/pis[j], after = j-1))
    })
    times = times[rows,]
    times = matrix(times, nrow = length(rows))
    rownames(times) = rows; colnames(times) = columns
    return(times)
  }else{ #User input transition rate matrix
    Q = args$transrate
    n = nrow(Q)
    pis = stationary(type = 'continuous', transrate = Q)$stationary.distribution
    times = sapply(columns, function(j){
      tmp = Q[,-j]; tmp = tmp[-j,]
      times = solve(-tmp, rep(1,n-1))
      return(append(times, 1/pis[j], after = j-1))
    })
    times = times[rows,]
    times = matrix(times, nrow = length(rows))
    rownames(times) = rows; colnames(times) = columns
    return(times)
  }
}

#Find the expected time for infinite states
timeinf = function(rows, columns, tol = 1e-6, ...){
  args = list(...)
  if(min(rows) < 1 || min(columns) < 1) stop('Start states at 1')
  
  truncate = function(K, transrate){
    Q = sapply(1:K, function(i){
      sapply(1:K, function(j){
        transrate(i,j)
      })
    })
    Q = t(Q)
    return(Q)
  }
  k = max(rows, columns)
  if(hasArg(transrate)){ #User input transition rates
    Q = truncate(k, args$transrate)
    
    #Check if the chain will stop or not
    check = stationary(type = 'continuous', transrate = Q)$stationary.distribution
    if(check[length(check)] > check[length(check)-1] && any(abs(diff(check)) >= 1e-10)) {
      if(min(rows) >= min(columns)) stop('no finite steps if i <= j')
      else warning('no finite steps if i <= j')
    }
    if(all(abs(diff(check)) < 1e-10)){
      if(any(rows %in% columns)) stop('no finite steps for i = i')
      else if(any(sapply(rows, function(i) i > columns))) stop('no finite steps for i < j')
      else warning('no finite steps for i = i')
    } 
    #Pass checking
    expecttime.old = timefin(rows = rows, columns = columns, transrate = Q)
    Q = truncate(k+50, args$transrate)
    expecttime.new = timefin(rows = rows, columns = columns, transrate = Q)
    error = norm(expecttime.old - expecttime.new, type = 'O')
    k = k+50
    
    while(error > tol){
      expecttime.old = expecttime.new
      k = k + 50
      Q = truncate(k, args$transrate)
      expecttime.new = timefin(rows=rows, columns=columns, transrate = Q)
      error = norm(expecttime.old - expecttime.new, type = 'O')
    }
    return(expecttime.new)
  }else{ #User input functions for probability matrix and holding rates
    rate.truncate = function(K, rates){
      qi = sapply(1:K, function(i){
        rates(i)
      })
      return(qi)
    }
    pij = truncate(k, args$pijdef)
    qi = rate.truncate(k, args$qidef)
    
    #Check if the chain will stop or not
    check = stationary(pijdef = pij, type = 'continuous', qidef = qi)$stationary.distribution
    if(check[length(check)] > check[length(check)-1] && any(abs(diff(check)) >= 1e-10)) {
      if(min(rows) >= min(columns)) stop('no finite steps if i <= j')
      else warning('no finite steps if i <= j')
    }
    if(all(abs(diff(check)) < 1e-10)){
      if(any(rows %in% columns)) stop('no finite steps for i = i')
      else if(any(sapply(rows, function(i) i > columns))) stop('no finite steps for i < j')
      else warning('no finite steps for i = i')
    } 
    #Pass checking
    expecttime.old = timefin(pijdef = pij, rows = rows, columns = columns, qidef = qi)
    pij = truncate(k+50, args$pijdef)
    qi = rate.truncate(k+50, args$qidef)
    expecttime.new = timefin(pijdef = pij, rows = rows, columns = columns, qidef = qi)
    error = norm(expecttime.old - expecttime.new, type = 'O')
    k = k+50
    
    while(error > tol){
      expecttime.old = expecttime.new
      k = k + 50
      pij = truncate(k, args$pijdef)
      qi = rate.truncate(k, args$qidef)
      expecttime.new = timefin(pijdef = pij, rows = rows, columns = columns, qidef = qi)
      error = norm(expecttime.old - expecttime.new, type = 'O')
    }
    return(expecttime.new)
  }
}