\name{umpbt.oneProp}
\alias{umpbt.oneProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT alternative used in the MSPRT in one-sample proportion tests
}
\description{
This function finds the UMPBT alternative, a 2-points mixture distribution, used in a MSPRT design in one-sample proportion tests (Table 1 in the main article). This is obtained by matching the rejection region of the UMPBT to that of the randomized fixed design test. For more details please refer to the following references.

This is a slight modification of the UMPBT point alternative as originally defined in Johnson (2013). This point alternative can be calculated using \code{point.umpbt.oneProp()}.
}
\usage{
umpbt.oneProp(side = "right", type1 = 0.005, n, null = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{type1}{
a numeric in \code{(0,1)}; the prespecified Type 1 error probability.

\strong{Default:} \code{0.005}.
}
  \item{n}{
a positive integer; number of samples to be used.
}
  \item{null}{
a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple null hypothesis.

\strong{Default} is \code{0.5}.
}
}
\details{
The UMPBT alternative used in the MSPRT in one-sample proportion tests is a mixture distribution with two points. This function returns those two points and the mixing probability.
}
\value{
Returns a list with the following two components:

\item{u}{a numeric vector of length 2; these are the two points of the mixture distribution.}

\item{psi}{a numeric in \code{(0,1)}; mixing probability corresponding to the first component of \code{u}.}
}
\references{
MSPRT: main article and supplemental information

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741

Johnson, Valen E., Revised standards for statistical evidence., Proceedings of the National Academy of Sciences, 16, 1945. (Specifically, it's supplemental file)

}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
# the point UMPBT alternative according to Johnson (2013)
point.umpbt.oneProp(n= 60, null= .2)
point.umpbt.oneProp(side= "left", n= 60, null= .2)

# the UMPBT alternative used by the MSPRT in these cases
umpbt.oneProp(n= 60, null= .2)
umpbt.oneProp(side= "left", n= 60, null= .2)
}

