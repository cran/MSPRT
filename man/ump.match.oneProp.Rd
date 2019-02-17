\name{ump.match.oneProp}
\alias{ump.match.oneProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding the "evidence threshold (\eqn{\delta})" in one-sample proportion tests in a fixed design
}
\description{
In one-sample proportion tests in a fixed design, this function solves for \eqn{\delta} by matching the rejection region from the UMPBT with that of the corresponding UMP (fixed design) test. Basically, this solves equation (20) in the supplemental information.
}
\usage{
ump.match.oneProp(side = "right", type1 = 0.005, n, p0 = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{type1}{
a numeric in \code{(0,1)}; prespecified Type 1 error probability.

\strong{Default:} \code{0.005}.
}
  \item{n}{
a positive numeric (integer); sample size to be used.
}
  \item{p0}{
a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple null hypothesis.

\strong{Default} is \code{0.5}.
}
}

\references{
MSPRT: Supplemental information

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
ump.match.oneProp(n= 60, p0= .2)
}

