\name{obj.func.ber}
\alias{obj.func.ber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Objective function for determining the UMPBT point alternative in a proportion test
}
\description{
The \eqn{h(p, \delta)} function as in the supplemental file. Given a \eqn{\delta}, we get the UMPBT(\eqn{\delta}) alternative by optimizing this function.
}
\usage{
obj.func.ber(p, delta, n.obs, p0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{a numeric in \code{(0,1)}; denotes a value of the binomial proportion}
  \item{delta}{a positive numeric; corresponding to this \eqn{\delta}, UMPBT(\eqn{\delta}) alternative is obtained}
  \item{n.obs}{a positive numeric (integer); number of samples to be used}
  \item{p0}{a numeric in \code{(0,1)}; the hypothesized value of the proportion under the simple null}
}

\value{Returns a numeric which is the value of the objective function}
\references{
MSPRT: Supplemental information

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{obj.func.ber(p= .5, delta= 25, n.obs= 60, p0= 0.2) }
