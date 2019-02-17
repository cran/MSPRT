\name{error.summary}
\alias{error.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error probabilities of a MSPRT
}
\description{
This function post processes the outputs of \code{overshoot.oneProp()}, \code{overshoot.oneZ()}, \code{overshoot.oneT()}, \code{overshoot.twoZ()} and \code{overshoot.twoT()}, and then calculates the Type 1 or Type 2 error probabilities of a MSPRT.
}
\usage{
error.summary(error.type, delta, root, count, inconclusive.vec, R,
              type1, type2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{error.type}{
a character;

            ``\code{type1}'' for computing Type 1 error probability.
            
            ``\code{type2}'' for computing Type 2 error probability.
}
  \item{delta}{
a positive numeric; the termination threshold.
}
  \item{root}{
a numeric; given this, (error probability - \code{root}) is returned.
}
  \item{count}{
a numeric (integer);

same as the ``\code{count}'' output from \code{overshoot.oneProp()}, \code{overshoot.oneZ()}, \code{overshoot.oneT()}, \code{overshoot.twoZ()} and \code{overshoot.twoT()} of the same \code{error.type}.
}
  \item{inconclusive.vec}{
a numeric vector;

same as the ``\code{inconclusive.vec}'' output from \code{overshoot.oneProp()}, \code{overshoot.oneZ()}, \code{overshoot.oneT()}, \code{overshoot.twoZ()} and \code{overshoot.twoT()} of the same \code{error.type}.
}
  \item{R}{
a numeric (integer); number of replications in Monte Carlo method that is used in \code{overshoot.oneProp()}, \code{overshoot.oneZ()}, \code{overshoot.oneT()}, \code{overshoot.twoZ()} and \code{overshoot.twoT()} of the same \code{error.type}.
}
  \item{type1}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 1 error of the MSPRT.
}
  \item{type2}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 2 error of the MSPRT.
}
}

\details{
Suppose we have the Monte Carlo outputs from \code{overshoot.oneProp()} or \code{overshoot.oneZ()} or \code{overshoot.oneT()} or \code{overshoot.twoZ()} or \code{overshoot.twoT()} corresponding to some \code{error.type}. Given these, the function returns the difference (error probability - \code{root}) for the MSPRT whose termination threshold is \code{delta}.
}
\value{
a numeric; the difference (error probability of the MSPRT - \code{root}).
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}
