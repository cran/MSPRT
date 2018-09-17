\name{point.umpbt.ber}
\alias{point.umpbt.ber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT point alternative in case of a proportion test
}
\description{
This function finds the UMPBT point alternative in case of a proportion test as originally defined in Johnson (2013). This is obtained by matching the rejection region of the UMPBT to that of the UMP (or fixed design) test.

\code{find.umpbt.ber()} is a slight modification of this.
}
\usage{
point.umpbt.ber( side = "right", type1 = 0.005, n.obs, null = 0.5)
}
%- maybe also 'usage' for other objects documented here.
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
point.umpbt.ber(n.obs= 60)
point.umpbt.ber(side= "left", n.obs= 60)
}

