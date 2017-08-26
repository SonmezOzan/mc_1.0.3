\name{plotstep}
\alias{plotstep}

\title{
  Plot heatmap of expected steps
}

\description{
Plots the expected steps for finite state irreducible Markov chain in a heatmap
}

\usage{
plotstep(pijdef, type, init = 1:nrow(pijdef), dest = 1:ncol(pijdef), ...)
}

\arguments{
  \item{pijdef}{The transition probabilities in matrix form. }
  \item{type}{Type of Markov chain, either 'discrete' or 'continuous'.}
  \item{init}{Set of desired initial states (default to all states).}
  \item{dest}{Set of desired destination states (default to all states).}
  \item{...}{Additional arguments for continuous Markov chain (see details).}
}


\details{ 
  
This function plots the heatmap of the expected steps of a finite state irreducible Markov chain. The user can supply the desired initial states and/or destination states using the arguments \emph{init} and/or \emph{dest} respectively. The function will plot the heatmap for the expected steps between states if these arguments are not supplied.

For continuous type markov chain, users must supply either both probability transition matrix (\emph{pijdef}) and holding rates (\emph{qidef}), or the transition rate matrix (\emph{transrate}) itself. See ?expectstep for more details.
}

\value{
A plot of the expected steps for the desired initial and destination states in the form of a heatmap.
}

\author{
Rex Cheung <rccheung@ucdavis.edu>, Teresa Filshtein <teresa.filshtein@gmail.com>, Norm Matloff <matloff@cs.ucdavis.edu>, and Ozan Sonmez <osonmez@ucdavis.edu>
}


\examples{
#Discrete Markov chain
P = matrix(rep(0.5, 9), ncol = 3)
P[1,3] = 0; P[2,2] = 0; P[3,1] = 0
plotstep(P, 'discrete')
plotstep(P, 'discrete', init = c(1,3))
plotstep(P, 'discrete', init = c(1,3), dest = c(2,3))

#Continuous Markov chain
#With probability transition matrix and holding rates
qidef = c(0.25, 0.175, 0.08)
pijdef = matrix(c(0, (1/20)/((1/20)+(1/8)), 0, 1, 0, 1, 0, (1/8)/((1/20)+(1/8)), 0), nrow = 3)
plotstep(pijdef, type = 'continuous', init = 1:2, dest = 1:3, qidef = qidef)
#With transition rates
Q = matrix(c(-0.25,0.05,0,0.25,-0.175,0.08,0,0.125,-0.08), nrow = 3)
plotstep(type = 'continuous', init = 1:2, dest = 1:3, transrate = Q)
}


