\name{absorb}
\alias{absorb}
\alias{absorb.mc}


\title{
Absorbing Markov Chain
}

\description{
Finds key features of a Discrete Finite absorbing Markov Chain. Future updates will include the infinite case}

\usage{
absorb.mc(pijdef, type, tol = 1e-06, ...)
}

\arguments{
   \item{pijdef}{The transition probabilities, either in matrix form or a function For now only matrix form, until the infinite case is incorporated }
   \item{type}{Type of markov chain, either 'discrete' or 'continuous'.}
   \item{tol}{A positive scalar for error tolerance for infinite markov chain approximation.}
   \item{...}{Additional argument for continuous markov chain type.}
}


\details{ 

This function generates key features and common output for an absorbing Markov Chain. An error is returned if the input matrix is not absorbing. A Markov Chain is an absorbing chain if 1) At least one state is absorbing 2) All non-absorbing states are transient.
}
\value{
Object of class "mc", with components
\item{pijdef}{The original transition matrix}
\item{absorb}{TRUE,Confirms it is an absorbing chain}
\item{astates}{absorbing states}
\item{tstates}{transient states}
\item{canonicalForm}{canonical form of the transition matrix}
\item{steady.state}{stationary distribution}
\item{FundamentalMatrix}{fundamental matrix}
\item{Qmat}{Q matrix, the sub matrix of transient to transient states}
\item{Rmat}{R matrix, the sub matrix of transient to absorbing states}
\item{ExpectedHits}{Expected number of hits in state i before getting absorbed}
\item{AbsorbProb}{A matrix of absorption probalities. Columns represent the absorbing states, rows represent transient states. The i,j entry represents the probablilty of being absorbed into absorbing state j given you are in state i}
}

\author{
Teresa Filshtein <teresa.filshtein@gmail.com>,Ozan Sonmez <osonmez@ucdavis.edu>, Rex Cheung <rccheung@ucdavis.edu>, and Norm Matloff <matloff@cs.ucdavis.edu> 
}


\examples{
#Discrete Chain
#Finite state


drunkard = t(matrix(c(1,rep(0,4),.5,0,.5,0,0,0,.5,0,.5,0,0,0,.5,0,.5,
                    rep(0,4),1),nrow = 5))
absorb.mc(drunkard)$ExpectedHits
absorb.mc(drunkard)$AbsorbProb
}




