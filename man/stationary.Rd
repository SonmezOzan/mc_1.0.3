\name{stationary}
\alias{stationary}
\alias{findpil.fin}
\alias{findpil.inf}
\alias{findpicont.fin}
\alias{findpicont.inf}

\title{
Stationary Distribution
}

\description{
Finds stationary distribution of an irreducible Markov chain. 
}

\usage{
stationary(pijdef, type, tol = 1e-06, ...)
findpil.fin(pijdef)
findpil.inf(pijdef, tol = 1e-06)
findpicont.fin(...)
findpicont.inf(tol = 1e-6,...)
}

\arguments{
   \item{pijdef}{The transition probabilities, either in matrix form (for finite states Markov chain) or a function (for infinite states Markov chain). }
   \item{type}{Type of Markov chain, either 'discrete' or 'continuous'.}
   \item{tol}{A positive scalar for error tolerance for infinite state Markov chain approximation.}
   \item{...}{Additional argument for continuous type Markov chain (see details).}
}


\details{ 

This function finds the stationary distribution of a given transition probability definition. User should input a matrix defining the transition probabilities for a finite state Markov chain, or a function specifying the probabilities for infinite state Markov chain (see example). 

The stationary distribution for the infinite states Markov chain is calculated via approximation. Each time we assume a finite number of states and calculate the stationary distributions, then increase the number of states and recalculate the stationary distributions, until the difference in between the distributions is less than the tolerate level.

For infinite state Markov chain, the forward probability should be smaller than the backward probability (except at the initial state), i.e. p_{ij} < p_{ik} for j > k, i > 1. This is to guarantee the chain will converge. An error will return if this condition is violated.

If the chain is of continuous type, user has the option of either supplying the holding rate at each state using the argument \emph{qidef}, or simply supplying the transition rate matrix using the argument \emph{transrate}. If supplying holding rates, the probability transition matrix must also be given (see example below).
}

\value{
Object of class "mc", with components 
\item{pijdef}{The input transition probability definition}
\item{stationary.distribution}{A numeric vector of probabilities for the stationary distribution}
\item{holding.rates}{If the chain is of continuous type and user supplied holding rates, it will also return the holding rates}
\item{transrate}{If the chain is of continuous type and user supplied transition rates, it will also return the transition rates}
}

\author{
Rex Cheung <rccheung@ucdavis.edu>, Teresa Filshtein <teresa.filshtein@gmail.com>, Norm Matloff <matloff@cs.ucdavis.edu>, and Ozan Sonmez <osonmez@ucdavis.edu>
}


\examples{
##Discrete Markov Chains
#Finite state
P = matrix(rep(0.5, 9), ncol = 3)
P[1,3] = 0; P[2,2] = 0; P[3,1] = 0
stationary(P, type = 'discrete')

#Infinite states
pijdef <- function(i,j){
 if (i == 1 && j == 2) return(0.3)
 if (i == 1 && j == 1) return(0.7)
 if ((i-j) == -1) return(0.3)
 if ((i-j) == 1) return(0.7)
 0
}
stationary(pijdef, type = 'discrete')

\dontrun{
##Continuous Markov Chains
#Finite States
#With probability transition matrix and holding rates
qidef = c(0.25, 0.175, 0.08)
pijdef = matrix(c(0, (1/20)/((1/20)+(1/8)), 0, 1, 0, 1, 0, (1/8)/((1/20)+(1/8)), 0), nrow = 3)
stationary(pijdef = pijdef, type = 'continuous', qidef = qidef)
#With transition rates
Q = matrix(c(-0.25, 0.25, 0, 0.05, -0.175, .125, 0, 0.08, -0.08), nrow = 3)
stationary(type = 'continuous', transrate = Q)

#Infinite States (Birth death process)
#With transition rates
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
#The diagonal elements are the negative of holding rates at each state
stationary(type = 'continuous', transrate = transrate, tol = 1e-1)

#With probability transition matrix and holding rates
qidef = function(i){
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
stationary(type = 'continuous', pijdef = pijdef, qidef = qidef, tol = 1e-1)
}
}


