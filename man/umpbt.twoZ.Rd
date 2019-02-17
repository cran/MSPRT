\name{umpbt.twoZ}
\alias{umpbt.twoZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT alternative in two-sample Z-tests
}
\description{
This function finds the UMPBT alternative in two-sample Z-tests. The value of hypothesized parameter, the difference between the population means of Group-2 and Group-1, under the simple null hypothesis is set at 0.
}
\usage{
umpbt.twoZ(side = "right", type1 = 0.005, n1, n2, sigma0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{type1}{
a numeric in \code{(0,1)}; prespecified Type 1 error probability.

\strong{Default:} \code{0.005}.
}
  \item{n1}{
a positive numeric (integer); sample size from Group-1 to be used.
}
  \item{n2}{
a positive numeric (integer); sample size from Group-2 to be used.
}
  \item{sigma0}{
a positive numeric; the known common population standard deviation.

\strong{Default:} 1.
}
}

\value{
Returns a numeric which is the UMPBT alternative in the two-sample Z-test.
}

\references{
MSPRT: main article and supplemental information

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741

Johnson, Valen E., Revised standards for statistical evidence., Proceedings of the National Academy of Sciences, 16, 1945. (Specially it's supplemental file)
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
umpbt.twoZ(n1= 60, n2= 60)
}

