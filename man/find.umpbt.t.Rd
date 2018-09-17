\name{find.umpbt.t}
\alias{find.umpbt.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT alternative in a T-test
}
\description{
This function finds the approximate data dependent UMPBT point alternative (Table 1 of the main article) in case of a T-test. For more details please refer to the supplemental information.
}
\usage{
find.umpbt.t(side = "right", type1 = 0.005, n.obs, null = 0, obs, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{
a character;

direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"};

\strong{Default} is \code{"right"}
}
  \item{type1}{
a numeric in \code{(0,1)}; the prespecified Type 1 error;

\strong{Default} is \code{0.005}
}
  \item{n.obs}{
a positive integer; number of samples to be used
}
  \item{null}{
a numeric; hypothesized value of population mean under the simple null hypothesis

\strong{Default} is \code{0}
}
  \item{obs}{
a numeric vector; a vector of the observations based on which the alternative is calculated;

If the alternative is required at step 5 when we are implementing a MSPRT, we need to provide all the data observed until that step in the order they were observed;

can be missing if \code{s} is provided.
}
  \item{s}{
a positive numeric; the sample standard deviation (sd) (of divisor \code{(n-1)}) of the \code{obs};

can be missing if \code{obs} is provided.
}
}

\details{
We need at least one of \code{obs} or \code{s}.
}

\value{
Returns a real numeric which is the UMPBT alternative in the T-test.
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
find.umpbt.t(n.obs= 60, obs= x.seq)

## providing the sd of x.seq
find.umpbt.t(n.obs= 60, s= sd(x.seq) )
}

