\name{umpbt.twoT}
\alias{umpbt.twoT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT alternative in two-sample T-tests
}
\description{
This function finds the approximate data dependent UMPBT alternative in two-sample T-tests. The value of hypothesized parameter, the difference between the population means of Group-2 and Group-1, under the simple null hypothesis is set at 0.
}
\usage{
umpbt.twoT(side = "right", type1 = 0.005, n1, n2, obs1, obs2, s)
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
  \item{obs1}{
a numeric vector; the vector of observations from Group-1 based on which the alternative is calculated.
}
  \item{obs2}{
a numeric vector; the vector of observations from Group-2 based on which the alternative is calculated.
}
  \item{s}{
a positive numeric; the pooled sample standard deviation (sd) based on \code{obs1} and \code{obs2}.

Can be missing if \code{obs1} and \code{obs2} are provided.
}
}

\details{
We need either \code{obs1} & \code{obs2}, or \code{s}. 

When we are implementing the MSPRT in a two-sample T-test and we need the alternative at Step-5, this function requires all the data from both groups observed until that step in the order they were observed.
}

\value{
Returns a numeric which is the UMPBT alternative in the two-sample T-test.
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
# suppose we get data from both groups after each observation
# a simulated ordered data at Step-30
x1.seq = rnorm(30,2,1.5)
x2.seq = rnorm(30,2,1.5)

## UMPBT alternative at step-30
# providing the data x1.seq & x2.seq
umpbt.twoT(n1= 60, n2= 60, obs1= x1.seq, obs2= x2.seq)
}

