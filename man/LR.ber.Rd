\name{LR.ber}
\alias{LR.ber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood ratio in a test for binomial proportion
}
\description{
Given a simple null and a simple alternative hypotheses, this function calculates the likelihood ratio (LR) in favor of the alternative based on an observed data in case of a test for binomial proportion.
}

\usage{
LR.ber(m, suff.stat, null = 0.5, alt)
}

\arguments{
  \item{m}{
a postive numeric (integer); number of samples for computing the LR
}
  \item{suff.stat}{
a postive numeric (integer); value of the sufficient statistic based on \code{m} observed data;

in this case the sufficient statistic is the total no. of successes out of \code{m} observations
}
  \item{null}{
a numeric in \code{(0,1)}; hypothesized value of the binomial proportion under the simple null

\strong{Default} is 0.5.
}
  \item{alt}{
a numeric in \code{(0,1)}; hypothesized value of the binomial proportion under the simple alternative
}
}

\value{
Returns a numeric denoting the LR in favor of the \code{alt} in the test for binomial proportion based on \code{m} observations.
}

\references{
MSPRT: supplemental information
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
LR.ber(m= 60, suff.stat= 48, null= 0.2, alt= 0.5)
}
