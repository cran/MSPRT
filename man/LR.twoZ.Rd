\name{LR.twoZ}
\alias{LR.twoZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood ratio in two-sample Z-tests
}
\description{
Given the simple null equals to 0 and a simple alternative hypotheses, this function calculates the likelihood ratio (LR) in favor of the alternative based on an observed data in two-sample Z-tests.
}
\usage{
LR.twoZ(m1, m2, suff.stat1, suff.stat2, alt, sigma0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{
a postive numeric (integer); number of samples from Group-1 for computing the LR.
}
  \item{m2}{
a postive numeric (integer); number of samples from Group-2 for computing the LR.
}
  \item{suff.stat1}{
numeric; the value of sufficient statistic based on \code{m1} samples from Group-1.

In this case, the sufficient statistic is the mean of those \code{m1} observations.
}
  \item{suff.stat2}{
numeric; the value of sufficient statistic based on \code{m2} samples from Group-2.

In this case, the sufficient statistic is the mean of those \code{m2} observations.
}
  \item{alt}{
a numeric; hypothesized value under the simple alternative
}
  \item{sigma0}{
a positive numeric; the known common population standard deviation (sd)

\strong{Default} is 1.
}
}


\value{
Returns a numeric denoting the LR in favor of \code{alt} in the two-sample Z-test based on the observations.
}

\references{
MSPRT: supplemental information
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
LR.twoZ(m1= 60, m2= 50, suff.stat1= 20, suff.stat2= 19, alt= 1, sigma0= 1.2)
}
