\name{error.summary}
\alias{error.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error probability of a MSPRT
}
\description{
After any of the functions \code{overshoot.ber()}, \code{overshoot.norm()} or \code{overshoot.t()} is executed, this function post processes it's outputs and then calculates the Type 1 (or Type 2) error probability of a MSPRT given a value of the termination threshold.
}
\usage{
error.summary(error.type, delta, root, count, inconclusive.vec, R,
              type1, type2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{error.type}{
a character;

            ``\code{type1}'' for computing Type 1 error probability
            
            ``\code{type2}'' for computing Type 2 error probability
}
  \item{delta}{
a positive numeric; a possible value of the termination threshold
}
  \item{root}{
a numeric; based on this (error probability - \code{root}) is returned
}
  \item{count}{
a numeric (integer);

same as the ``\code{count}'' output from \code{overshoot.ber()}, \code{overshoot.norm()} or \code{overshoot.t()} with the same \code{error.type}
}
  \item{inconclusive.vec}{
a numeric vector;

same as the ``\code{inconclusive.vec}'' output from \code{overshoot.ber()}, \code{overshoot.norm()} or \code{overshoot.t()} with the same \code{error.type}
}
  \item{R}{
a numeric (integer); number of replications in Monte Carlo method that is used in \code{overshoot.ber()}, \code{overshoot.norm()} or \code{overshoot.t()} with the same \code{error.type}
}
  \item{type1}{
a numeric in \code{(0,1)}; the prespecified Type 1 error probability used in \code{overshoot.ber()}, \code{overshoot.norm()} or \code{overshoot.t()}
}
  \item{type2}{
a numeric in \code{(0,1)}; the prespecified Type 2 error probability used in \code{overshoot.ber()}, \code{overshoot.norm()} or \code{overshoot.t()}
}
}

\details{
Suppose we have the Monte Carlo method outputs from \code{overshoot.ber()}, \code{overshoot.norm()} or \code{overshoot.t()} corresponding to some \code{error.type}. Given this, the function calculates the difference (error probability - \code{root}) of a MSPRT where \code{delta} is used as a value of the termination threshold.
}
\value{
a numeric; the difference (estimated error probability - \code{root})
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}
