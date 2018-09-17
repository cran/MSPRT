\name{find.umpbt.ber}
\alias{find.umpbt.ber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT alternative used in the MSPRT in case of a proportion test
}
\description{
This function finds the UMPBT alternative, the 2-points mixture distribution, used by the MSPRT design in case of a proportion test. This is obtained by matching the rejection region of the UMPBT to that of the randomized fixed design test. For more details please refer to the supplemental information.

This is a slight modification of the UMPBT point alternative as originally defined in Johnson (2013). The original point alternative can be calculated by using \code{point.umpbt.ber()}.
}
\usage{
find.umpbt.ber(side = "right", type1 = 0.005, n.obs, null = 0.5)
}

\arguments{
  \item{side}{
a character;

direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"};

\strong{Default} is \code{"right"}
}
  \item{type1}{
a numeric in \code{(0,1)}; the prespecified Type 1 error

\strong{Default} is \code{0.005}
}
  \item{n.obs}{
a positive integer; number of samples to be used
}
  \item{null}{
a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple null hypothesis

\strong{Default} is \code{0.5}
}
}
\details{
The UMPBT alternative used in the MSPRT is a mixture distribution with two points. This function returns those two points and the mixing probability.
}
\value{
Returns a list with the following two components:

\item{u}{a numeric vector of length 2; these are the two points of the mixture distribution}

\item{psi}{a numeric in \code{(0,1)}; mixing probability corresponding to the first component of \code{u}}
}
\references{
MSPRT: main article and supplemental information

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741

Johnson, Valen E., Revised standards for statistical evidence., Proceedings of the National Academy of Sciences, 16, 1945.(Specially it's supplemental file)

}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
find.umpbt.ber(side= "right", n.obs= 60, null= .2)
find.umpbt.ber(side= "left", n.obs= 60, null= .2)
}
