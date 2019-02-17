\name{LR.oneProp}
\alias{LR.oneProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood ratio in one-sample proportion tests
}
\description{
Given a simple null and a simple alternative hypotheses, this function calculates the likelihood ratio (LR) in favor of the alternative based on an observed data in one-sample proportion tests.
}
\usage{
LR.oneProp(m, suff.stat, null = 0.5, alt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
a postive numeric (integer); number of samples for computing the LR.
}
  \item{suff.stat}{
a postive numeric (integer); the value of sufficient statistic based on \code{m} observed data.

In this case, the sufficient statistic is the total no. of successes out of \code{m} observations.
}
  \item{null}{
a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple null.

\strong{Default:} 0.5.
}
  \item{alt}{
a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple alternative.
}
}

\value{
Returns a numeric denoting the LR in favor of \code{alt} in the one-sample proportion test based on \code{m} observations.
}

\references{
MSPRT: supplemental information
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
LR.oneProp(m= 60, suff.stat= 48, null= 0.2, alt= 0.5)
}
