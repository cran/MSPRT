\name{ump.match.ber}
\alias{ump.match.ber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding the "evidence threshold (\eqn{\delta})" in a proportion test
}
\description{
This function solves for \eqn{\delta} by matching the rejection region from the UMPBT with that of the corresponding fixed design test. Basically this solves equation (20) of the supplemental information.
}
\usage{
ump.match.ber( side = "right", type1 = 0.005, n.obs, p0 = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{
a character;

direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"};

\strong{Default} is \code{"right"}.
}
  \item{type1}{
a numeric in \code{(0,1)}; the prespecified Type 1 error;

\strong{Default} is \code{0.005}.
}
  \item{n.obs}{
a positive integer; number of samples to be used
}
  \item{p0}{
a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple null hypothesis;

\strong{Default} is \code{0.5}.
}
}


\value{
Returns a numeric which is the value of \eqn{\delta}, the "evidence threshold", as defined in Johnson (2013). The value of \eqn{\delta} is such that the rejection region from the UMPBT matches with that of the corresponding fixed design test.
}
\references{
MSPRT: Supplemental information

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
ump.match.ber(n.obs= 60, p0= .2)
}
