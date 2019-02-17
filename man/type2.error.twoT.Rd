\name{type2.error.twoT}
\alias{type2.error.twoT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Type 2 error function of two-sample T-tests in a fixed design
}
\description{
For a two-sample T-test in a fixed design, this function evaluates (Type 2 error - a constant) at a specified alternative value of the hypothesized parameter (the difference between the population means of Group-2 and Group-1). The hypothesized value under the simple null hypothesis is set at 0.
}
\usage{
type2.error.twoT(alt, side = "right", n1, n2, type1 = 0.005, root = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alt}{
a numeric; the value of hypothesized parameter (consistent with \code{side}) where the Type 2 error of the fixed design test needs to be evaluated at.
}
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{n1}{
a positive numeric (integer); sample size from Group-1 to be used.
}
  \item{n2}{
a positive numeric (integer); sample size from Group-2 to be used.
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
In case of two-sample T-tests in a fixed design, this function evaluates the Type 2 error at a specified value \code{alt}. We can also use this function to do the other way round by exploiting the argument \code{root}; that is, given a Type 2 error \eqn{\beta}, we can obtain the \code{alt} where Type 2 error equals \eqn{\beta}. To do that, we need to substitute \code{root}=\eqn{\beta} in the argument of this function, and then solve it for \code{alt}. The function \code{find.alt()} in this package excatly does this.
}

\value{
If \code{root=k}, a numeric \code{(Type 2 error - k)} evaluated at the speified alternative value \code{alt} is returned.
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
## Type 2 error at an alternative value
type2.error.twoT(alt= 1.2, n1= 60, n2= 60)

## (Type 2 error - 0.5) at the same alternative value
type2.error.twoT(alt= 1.2, n1= 60, n2= 60, root = 0.5)
}

