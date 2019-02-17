\name{LR.twoT}
\alias{LR.twoT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayes factor in two-sample T-tests
}
\description{
Given the simple null equals to 0 and a simple alternative hypothesis, this function calculates the bayes factor (BF) in favor of the alternative based on an observed data in two-sample T-tests.
}
\usage{
LR.twoT(m1, m2, suff.stat1, suff.stat2, alt, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{
a postive numeric (integer); number of samples from Group-1 for computing the BF.
}
  \item{m2}{
a postive numeric (integer); number of samples from Group-2 for computing the BF.
}
  \item{suff.stat1}{
numeric; the value of sufficient statistic based on \code{m1} samples from Group-1.

In this case, the sufficient statistic is the mean of these \code{m1} observations.
}
  \item{suff.stat2}{
numeric; the value of sufficient statistic based on \code{m2} samples from Group-2.

In this case, the sufficient statistic is the mean of these \code{m2} observations.
}
  \item{alt}{
a numeric; the hypothesized value under the simple alternative.
}
  \item{s}{
a positive numeric; pooled sample standard deviation using \code{m1} observations from Group-1 and \code{m2} observations from Group-2.
}
}

\value{
Returns a numeric denoting the BF in favor of \code{alt} in the two-sample T-test based on the observations.
}

\references{
MSPRT: supplemental information
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
LR.twoT(m1= 60, m2= 50, suff.stat1= 20, suff.stat2= 19, alt= 1, s= 1)
}



