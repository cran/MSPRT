\name{type2.error.t}
\alias{type2.error.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Type 2 error function of the T-test in a fixed design
}
\description{
For a T-test in a fixed design, this evaluates (Type 2 error - a constant) at a specified value of the population mean. This value should be consistent with the direction of the alternative hypothesis.
}
\usage{
type2.error.t( alt, side = "right", null.val = 0, N, type1 = 0.005, root = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{
a character;

direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"};

\strong{Default} is \code{"right"}.
}
  \item{alt}{
a numeric; the value of the population mean (consistent with \code{side}) where the Type 2 error of the fixed design test will be evaluated
}
  \item{null.val}{
a numeric; hypothesized value of the population mean under a simple null hypothesis;

\strong{Default} is \code{0}.
}
  \item{N}{
a positive numeric (integer); sample size to be used
}
  \item{type1}{
a numeric in \code{(0,1)}; prespecified Type 1 error;

\strong{Default} is \code{0.005}.
}
  \item{root}{
a numeric; \strong{Default} is 0, in which case the Type 2 error is returned;

in general, \code{(}Type 2 error - \code{root)} is returned;

\strong{Default} is \code{0}.
}
}

\details{
In case of a T-test in a fixed design, this function evaluates the Type 2 error at a specified value \code{alt} in the composite alternative region. But we can also use this function to do the other way round by using the argument \code{root}; that is, given a Type 2 error \eqn{\beta}, we can obtain the value of the proportion in the composite alternative region which has Type 2 error \eqn{\beta}. For that we need to substitute root=\eqn{\beta} in the argument of this function, and then solve it for \code{alt}. The function \code{find.alt()} in this package excatly does this.
}
\value{
If \code{root=k}, this returns a numeric \code{(}Type 2 error - \code{k)} at the speified alternative value \code{alt}.
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
## Type 2 error at an alternative value
type2.error.t(alt= 1.2, N= 60)

## (Type 2 error - 0.5) at the same alternative value
type2.error.t(alt= 1.2, N= 60, root = 0.5)
}


