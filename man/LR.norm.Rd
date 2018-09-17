\name{LR.norm}
\alias{LR.norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood ratio in a Z-test
}
\description{
Given a simple null and a simple alternative hypotheses, this function calculates the likelihood ratio (LR) in favor of the alternative based on an observed data in case of a Z-test.}

\usage{
LR.norm(m, suff.stat, null = 0, alt, sigma0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
a postive numeric (integer); number of samples for computing the LR
}
  \item{suff.stat}{
a postive numeric (integer); value of the sufficient statistic based on \code{m} observed data;

in this case the sufficient statistic is the sum of \code{m} observations
}
  \item{null}{
a numeric; hypothesized value of population mean under the simple null

\strong{Default} is 0.
}
  \item{alt}{
a numeric; hypothesized value of population mean under the simple alternative
}
  \item{sigma0}{
a positive numeric; known standard deviation (sd) of the data in a Z-test;

\strong{Default} is 1.
}
}

\value{
Returns a numeric denoting the LR in favor of the \code{alt} in case of the Z-test based on \code{m} observations.
}

\references{
MSPRT: supplemental information
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
LR.norm(m= 60, suff.stat= 10.5, alt= 0.5)
}

