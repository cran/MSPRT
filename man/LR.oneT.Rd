\name{LR.oneT}
\alias{LR.oneT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayes factor in one-sample T-tests
}
\description{
Given a simple null and a simple alternative hypotheses, this function calculates the bayes factor (BF) in favor of the alternative based on an observed data in one-sample T-tests.}
\usage{
LR.oneT(m, suff.stat, null = 0, alt, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
a postive numeric (integer); number of samples for computing the BF.
}
  \item{suff.stat}{
numeric; the value of sufficient statistic based on \code{m} observed data.

In this case, the sufficient statistic is the sum of \code{m} observations.
}
  \item{null}{
a numeric; the hypothesized value of population mean under the simple null.

\strong{Default:} 0.
}
  \item{alt}{
a numeric; the hypothesized value of population mean under the simple alternative.
}
  \item{s}{
a positive numeric; sample standard deviation (with divisor \code{(m-1)}) of the \code{m} observations.
}
}

\value{
Returns a numeric denoting the BF in favor of \code{alt} in one-sample T-tests based on \code{m} observations.
}

\references{
MSPRT: supplemental information
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
LR.oneT(m= 60, suff.stat= 20.2, alt= 1.5, s= 1.2)
}

