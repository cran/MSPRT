\name{LR.t}
\alias{LR.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayes factor in a T-test
}
\description{
Given a simple null and a simple alternative hypotheses, this function calculates the bayes factor (BF) in favor of the alternative based on an observed data in case of a T-test.}

\usage{
LR.t(m, suff.stat, null = 0, alt, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
a postive numeric (integer); number of samples for computing the BF
}
  \item{suff.stat}{
a postive numeric (integer); value of the sufficient statistic based on \code{m} observed data;

in this case the sufficient statistic is the sum of \code{m} observations.
}
  \item{null}{
a numeric; hypothesized value of population mean under the simple null

\strong{Default} is 0.
}
  \item{alt}{
a numeric; hypothesized value of population mean under the simple alternative
}
  \item{s}{
a positive numeric; sample standard deviation (with divisor \code{(n-1)}) of the \code{m} observations
}
}

\value{
Returns a numeric denoting the BF in favor of the specified \code{alt} in case of T-test based on \code{m} observations.
}

\references{
MSPRT: supplemental information
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
LR.t(m= 60, suff.stat= 20.2, alt= 1.5, s= 1.2)
}
