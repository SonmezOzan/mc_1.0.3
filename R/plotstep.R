#Input Variables:
#pijdef: a matrix capturing the transition probabilities
#type: specify whether the chain is discrete or continuous
#init: set of interested initial states 
#dest: set of interested destination states

#Note: This function is for finite states and irreducible chains only



plotstep = function(pijdef, type, init = 1:nrow(pijdef), dest = 1:ncol(pijdef), ...){
 
  args = list(...)
  init = sort(unique(init)); dest = sort(unique(dest))
  if(type == 'discrete'){
    step.matrix = expectstep(pijdef, type, init, dest)$steps
    step.matrix = step.matrix[nrow(step.matrix):1, , drop = FALSE]
  }else{
    if(hasArg(transrate)){
      step.matrix = expectstep(type = 'continuous', rows = init, columns = dest, transrate = args$transrate)$time
      step.matrix = step.matrix[nrow(step.matrix):1, , drop = FALSE]
    }else{
      if(!hasArg(qidef)) stop('Missing holding rates')
      step.matrix = expectstep(pijdef, type = 'continuous', rows = init, columns = dest, qidef = args$qidef)$time
      step.matrix = step.matrix[nrow(step.matrix):1, , drop = FALSE]
    }
  }
  grid = expand.grid(x = 1:nrow(step.matrix), y = 1:ncol(step.matrix))
  grid$z = as.vector(step.matrix)
  print(levelplot(z ~ y*x, grid, xlab = 'Destination States', ylab = 'Initial States', main = 'Expected Steps/Time',
                  col.regions = heat.colors,
                  panel = function(...){
                    arg = list(...)
                    panel.levelplot(...)
                    panel.text(arg$x, arg$y, round(arg$z,2))
                  },
                  scales = list(x = list(at = 1:ncol(step.matrix), labels = dest), 
                                y = list(at = 1:nrow(step.matrix), labels = rev(init))),
                  aspect="fill"))
}


















