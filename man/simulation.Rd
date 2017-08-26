\name{simulation}
\alias{mc.simulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Markov Chain Simulation
}

\description{
Simulates discrete markov chain with state space}

\usage{
mc.simulation(pijdef, type, N, initial.state,...)
}

\arguments{
   \item{pijdef}{The transition probabilities, either in matrix form or a function For now only matrix form, until the infinite case is incorporated }
   \item{type}{Type of markov chain, either 'discrete' or 'continuous'.}
   \item{N}{Number of steps of the markov chain being simulated}
   \item{initial.state}{Initial state that the chain will start from}
   \item{...}{Additional argument for continuous markov chain type.}
}


\details{ 

This function will simulate discrete markov chain for a given transition probability matrix, if the initial state is not specified then it will be randomly generated from the state space.
}
\value{
\item{pijdef}{The original transition matrix}
\item{simulated.states}{N Simulated states based on the transition matrix}
}

\author{
Teresa Filshtein <teresa.filshtein@gmail.com>, Ozan Sonmez <osonmez@ucdavis.edu>, Rex Cheung <rccheung@ucdavis.edu>, and Norm Matloff <matloff@cs.ucdavis.edu> 
}


\examples{
t.matrix = t(matrix(c(0.1,0.2,0.7,0.4,0.5,0.1,0,0.4,0.6),3,3))
mc.simulation(t.matrix, "discrete", 100, 1)$simulated.states
}


