\name{effective.N}
\alias{effective.N}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determining the 'effective maximum sample size' in MSPRT for one-sample proportion tests
}
\description{
Because of the discreteness issue in one-sample proportion tests, it is not always 'effective' (the sense is explained in the Details section) to use just any value as a maximum sample size for designing a MSPRT. Suppose, we have a desired value for the maximum sample size. Given this, the function finds a value, defined as the effective maximum sample size, which should be used as the maximum sample size for designing the MSPRT (as \code{N.max} in \code{design.MSPRT()} ).
}
\usage{
effective.N(N, side = "right", type1 = 0.005, null = 0.5, plot.it = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
a numeric; a desired value of the maximum sample size for designing the MSPRT.
}
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{type1}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 1 error of the MSPRT.

\strong{Default:} \code{0.005}.
}
  \item{null}{
a numeric; the hypothesized value of proportion under the null hypothesis.

\strong{Default:} 0.5.
}
  \item{plot.it}{
a logical; if TRUE, returns a plot; otherwise it doesn't.

\strong{Default:} TRUE.
}
}
\details{
Suppose we are provided with a simple null hypothesis, and Type 1 & 2 error probabilities. Because of the discreteness issue in one-sample proportion tests in a fixed design, the fixed design alternative will not always decrease even if we are increasing the sample size. So, we first shortlist only those values from 1 to \code{N} which results in strictly decreasing UMPBT point alternatives (as is originally defined in Johnson (2013)). The effective maximum sample size is then chosen to be the maximum among those shorlisted value. So, obviously, the effective maximum sample size is either \code{N} or smaller than that.
}

\value{
Returns a numeric, the effective maximum sample size.
}
\references{
Main article and supplemental file of MSPRT

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
effective.N(N = 30, null = .2)
}
