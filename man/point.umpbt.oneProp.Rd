\name{point.umpbt.oneProp}
\alias{point.umpbt.oneProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT point alternative in one-sample proportion tests
}
\description{
This function finds the UMPBT point alternative, as defined in Johnson (2013), in case of one-sample proportion tests. This is obtained by matching the rejection region of the UMPBT to that of the UMP (or fixed design) test.
}
\usage{
point.umpbt.oneProp(side = "right", type1 = 0.005, n, null = 0.5)
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

\value{
Returns a numeric which is the UMPBT point alternative as defined in Johnson (2013).
}
\references{
Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
point.umpbt.oneProp(n= 60, null= .2)
point.umpbt.oneProp(side= "left", n= 60, null= .2)
}


