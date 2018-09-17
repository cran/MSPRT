\name{effective.N}
\alias{effective.N}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effective maximum sample size for a MSPRT in a proportion test
}
\description{
Because of the discreteness issue in a proportion test, it may not be effective if we use just any value as a maximum sample size to design a MSPRT. Given a desired value for the maximum sample size, this function finds the maximum sample size which should be used for the design.
}
\usage{
effective.N(N.max, side = "right", type1 = 0.005, null.val = 0.5, plot.it = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N.max}{
a numeric; a desired value of the maximum sample size for designing a MSPRT
}
  \item{side}{
a character;

direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"}.

\strong{Default} is \code{"right"}.
}
  \item{type1}{
a numeric in \code{(0,1)};

The probability at which we want to control the Type 1 error of the MSPRT.

\strong{Default} is \code{0.005}.
}
  \item{null.val}{
a numeric; 

The hypothesized value of the proportion under the null hypothesis.

\strong{Default} is 0.5.
}
  \item{plot.it}{
a logical; if TRUE, returns a plot. Otherwise not.

\strong{Default} is TRUE.
}
}
\details{
First we shortlist only those values within \code{N.max} at which the UMPBT point alternative, as is originally defined in Johnson (2013), strictly decreases. The final value is chosen to be the maximum among those shorlisted value.
}

\value{
Returns a numeric. This is the `effective' maximum sample size which should be used to design a MSPRT.
}
\references{
Main article and supplemental file of MSPRT

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
effective.N(N.max = 30, null.val = .2)

}
