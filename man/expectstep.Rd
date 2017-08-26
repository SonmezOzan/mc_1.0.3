\name{expectstep}
\alias{expectstep}
\alias{steps.fin}
\alias{steps.inf}
\alias{timefin}
\alias{timeinf}

\title{
Expected Steps/Time
}

\description{
Finds the expected steps/time from one state to another. 
}

\usage{
expectstep(pijdef, type, rows, columns, tol = 1e-6, ...)
steps.fin(pijdef, type, rows, columns)
steps.inf(pijdef, type, rows, columns, tol = 1e-06)
timefin(rows, columns, ...)
timeinf(rows, columns, tol = 1e-6, ...)
}

\arguments{
   \item{pijdef}{The transition probabilities, either in matrix form (for finite states Markov chain) or a function (for infinite states Markov chain). }
   \item{type}{Type of Markov chain, either 'discrete' or 'continuous'.}
   \item{rows}{A numeric vector of initial states of interest.}
   \item{columns}{A numeric vector of destination states of interest.}
   \item{tol}{A positive scalar for error tolerance for infinite state Markov chain approximation.}
   \item{...}{Additional argument for continuous type Markov chain (see details).}
}


\details{ 

This function finds the expected steps it takes to reach from one state to another in an irreducible Markov chain. User should input a matrix defining the transition probabilities for a finite state Markov chain, or a function specifying the probabilities for infinite state Markov chain (see example). 

The expected steps for the infinite state Markov chain is calculated via approximation. Each time we assume a finite number of states and calculate the expected number of steps, then increase the number of states and recalculate the expected steps, until the one norm of the difference in the step matrices is less than the tolerate level.

For infinite state Markov chain, the forward probability should be smaller than the backward probability (except at the initial states), i.e. p_{ij} < p_{ik} for j > k. This is to guarantee the chain will converge. If such condition is violated, only the expected steps from states i to states j for all i < j can be calculated, or an error will be returned otherwise. 

For continuous type Markov chain, user has the option of either supplying the holding rate at each state using the argument \emph{qidef}, or simply supplying the transition rate matrix using the argument \emph{transrate}. If supplying holding rates, the probability transition matrix must also be given (see example below). Note \emph{diag(transrate) = -qidef}.
}

\value{
Object of class "mc", with components 
\item{pijdef}{The input transition probability definition if supplied.}
\item{steps}{(For discrete chains) The expected steps in a matrix form, with rows specifying the initial states and columns specifying the destination states.}
\item{time}{(For continuous chains) The expected time in a matrix form, with rows specifying the initial states and columns specifying the destination states.}
\item{holding.rates}{If the chain is of continuous type and user supplied holding rates, it will also return the holding rates}
\item{transrate}{If the chain is of continuous type and user supplied transition rates, it will also return the transition rates}
}

\author{
Rex Cheung <rccheung@ucdavis.edu>, Teresa Filshtein <teresa.filshtein@gmail.com>, Norm Matloff <matloff@cs.ucdavis.edu>, and Ozan Sonmez <osonmez@ucdavis.edu>
}


\examples{
##Discrete Markov Chain
#Finite state
P = matrix(rep(0.5, 9), ncol = 3)
P[1,3] = 0; P[2,2] = 0; P[3,1] = 0
expectstep(P, type = 'discrete', rows = 1:3, columns = 1:3)

#Infinite states
pijdef <- function(i,j){
  if (i == 1 && j == 2) return(0.4)
  if (i == 1 && j == 1) return(0.6)
  if ((i-j) == -1) return(0.4)
  if ((i-j) == 1) return(0.6)
  0
}
expectstep(pijdef, type = 'discrete', rows = c(2,4,5), columns = c(2,3,7,9))

\dontrun{
##Continuous Markov Chain
#Finite States
#With probability transition matrix and holding rates
qidef = c(0.25, 0.175, 0.08)
pijdef = matrix(c(0, (1/20)/((1/20)+(1/8)), 0, 1, 0, 1, 0, (1/8)/((1/20)+(1/8)), 0), nrow = 3)
expectstep(pijdef, 'continuous', rows = 1:3, columns = 1:3, qidef = qidef)
#With transition rates
Q = matrix(c(-0.25,0.05,0,0.25,-0.175,0.08,0,0.125,-0.08), nrow = 3)
expectstep(type = 'continuous', rows = 1:3, columns = 1:3, transrate = Q)

#Infinite States (Birth death process)
#With Transition rates
transrate = function(i,j){
  if(i == 1 && j == 1) return(-1)
  if(i == 1 && j == 2) return(1)
  if(i == 2 && j == 1) return(1)
  if(i == 2 && j == 2) return(-4/3)
  if(i == 2 && j == 3) return(1/3)
  if(j - i == 1) return(1/(i+1)) #Foward
  if(j - i == -1) return(1/(i-1)) #Backward
  if(j - i == 0) return(-1/(i-1)-1/(i+1)) #Itself
  0
}
expectstep(pijdef = NULL, type = 'continuous', rows = 1:3, columns = 1:3, transrate = transrate, tol = 0.5)

#With probability transition matrix and holding rates
holding = function(i){
  if(i == 1) return(1)
  return(1/(i-1) + 1/(i+1))
}
pijdef = function(i,j){
  if(i == 1 && j == 2) return(1)
  if(i == 1 && j == 1) return(0)
  if(j - i == 1) return((1/(i+1))/((1/(i-1)) + (1/(i+1))))
  if(j - i == 0) return(0)
  if(j - i == -1) return((1/(i-1))/((1/(i-1)) + (1/(i+1))))
  0
}
expectstep(pijdef = pijdef, type = 'continuous', rows = 1:3, columns = 1:3, qidef = holding, tol = 0.5)
}



}


