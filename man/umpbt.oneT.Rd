\name{umpbt.oneT}
\alias{umpbt.oneT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT alternative in one-sample T-tests
}
\description{
This function finds the approximate data dependent UMPBT alternative (Table 1 in the main article) in one-sample T-tests. For more details please refer to the following references.
}
\usage{
umpbt.oneT(side = "right", type1 = 0.005, n, null = 0, obs, s)
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
  \item{obs}{
a numeric vector; the vector of observations based on which the alternative needs to be calculated. This is ignored if \code{s} is provided.
}
  \item{s}{
a positive numeric; the sample standard deviation (sd) (of divisor \code{(n-1)}) of the \code{obs}.

Can be missing if \code{obs} is provided.
}
}

\details{
We need either \code{obs} or \code{s}. 

When we are implementing the MSPRT in one-sample T-test and we need the alternative at Step-5, this function requires all the data observed until Step-5 in the order they were observed.
}

\value{
Returns a numeric which is the UMPBT alternative in the one-sample T-test.
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
# a simulated ordered data at step-30
x.seq = rnorm(30,2,1.5)

# UMPBT alternative at step-30

## providing the data x.seq
umpbt.oneT(n= 60, obs= x.seq)

## providing the sd of x.seq
umpbt.oneT(n= 60, s= sd(x.seq))
}

