\name{type2.error.oneT}
\alias{type2.error.oneT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Type 2 error function of one-sample T-tests in a fixed design
}
\description{
For a one-sample T-test in a fixed design, this function evaluates (Type 2 error - a constant) at a specified alternative value of population mean.
}
\usage{
type2.error.oneT(alt, side = "right", null = 0, n, type1 = 0.005, root = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alt}{
a numeric; the value of population mean (consistent with \code{side}) where the Type 2 error of the fixed design test needs to be evaluated at.
}
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{null}{
a numeric; the hypothesized value of population mean under a simple null hypothesis.

\strong{Default:} \code{0}.
}
  \item{n}{
a positive numeric (integer); sample size to be used.
}
  \item{type1}{
a numeric in \code{(0,1)}; prespecified Type 1 error probability.

\strong{Default:} \code{0.005}.
}
  \item{root}{
a numeric; when this is 0, the Type 2 error is returned;

in general, \code{(Type 2 error - root)} is returned

\strong{Default:} \code{0}.
}
}

\details{
In case of one-sample T-tests in a fixed design, this function evaluates the Type 2 error at a specified value \code{alt}. We can also use this function to do the other way round by exploiting the argument \code{root}; that is, given a Type 2 error \eqn{\beta}, we can obtain the \code{alt} where Type 2 error equals \eqn{\beta}. To do that, we need to substitute \code{root}=\eqn{\beta} in the argument of this function, and then solve it for \code{alt}. The function \code{find.alt()} in this package excatly does this.
}
\value{
If \code{root=k}, a numeric \code{(Type 2 error - k)} evaluated at the speified alternative value \code{alt} is returned.
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
## Type 2 error at an alternative value
type2.error.oneT(alt= 1.2, n= 60)

## (Type 2 error - 0.5) at the same alternative value
type2.error.oneT(alt= 1.2, n= 60, root = 0.5)
}

