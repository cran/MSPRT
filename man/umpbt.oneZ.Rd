\name{umpbt.oneZ}
\alias{umpbt.oneZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT alternative in one-sample Z-tests
}
\description{
This function finds the UMPBT alternative (Table 1 in the main article) in one-sample Z-tests. This is obtained by matching the rejection region of the UMPBT to that of the fixed design test. For more details please refer to the following references.
}
\usage{
umpbt.oneZ(side = "right", type1 = 0.005, n, null = 0, sigma0 = 1)
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
  \item{n}{
a positive numeric (integer); sample size to be used.
}
  \item{null}{
a numeric; the hypothesized value of population mean under the simple null hypothesis.

\strong{Default:} \code{0}.
}
  \item{sigma0}{
a positive numeric; the known population standard deviation.

\strong{Default:} 1.
}
}

\value{
Returns a numeric which is the UMPBT point alternative in the one-sample Z-test.
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
umpbt.oneZ(n= 60)
}

