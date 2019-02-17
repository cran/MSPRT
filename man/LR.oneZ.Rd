\name{LR.oneZ}
\alias{LR.oneZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood ratio in one-sample Z-tests
}
\description{
Given a simple null and a simple alternative hypotheses, this function calculates the likelihood ratio (LR) in favor of the alternative based on an observed data in one-sample Z-tests.
}
\usage{
LR.oneZ(m, suff.stat, null = 0, alt, sigma0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
a postive numeric (integer); number of samples for computing the LR.
}
  \item{suff.stat}{
numeric; the value of sufficient statistic based on \code{m} observed data.

In this case, the sufficient statistic is the sum of \code{m} observations
}
  \item{null}{
a numeric; the hypothesized value of population mean under the simple null.

\strong{Default:} 0.
}
  \item{alt}{
a numeric; the hypothesized value of population mean under the simple alternative.
}
  \item{sigma0}{
a positive numeric; the known population standard deviation (sd).

\strong{Default:} 1.
}
}

\value{
Returns a numeric denoting the LR in favor of \code{alt} in the one-sample Z-test based on \code{m} observations.
}

\references{
MSPRT: supplemental information
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
LR.oneZ(m= 60, suff.stat= 10.5, alt= 0.5)
}
