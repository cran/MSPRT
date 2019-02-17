\name{objfunc.oneProp}
\alias{objfunc.oneProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Objective function for determining the UMPBT point alternative in one-sample proportion tests
}
\description{
This is the \eqn{h(p, \delta)} function as in the supplemental file. Given a \eqn{\delta}, we optimize this function to get the UMPBT(\eqn{\delta}) alternative for one-sample proportion tests in a fixed design.
}
\usage{
objfunc.oneProp(p, delta, n, p0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{a numeric in \code{(0,1)}; the value of proportion.}
  \item{delta}{a positive numeric; corresponding to this \eqn{\delta}, the UMPBT(\eqn{\delta}) alternative is obtained.}
  \item{n}{a positive numeric (integer); number of samples to be used.}
  \item{p0}{a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple null.}
}

\value{Returns a numeric which is the value of the objective function.}
\references{
MSPRT: Supplemental information

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
objfunc.oneProp(p= .5, delta= 25, n= 60, p0= 0.2)
}
