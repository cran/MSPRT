\name{find.umpbt.norm}
\alias{find.umpbt.norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The UMPBT alternative in a Z-test
}
\description{
This function finds the UMPBT point alternative (Table 1 of the main article) in case of a Z-test. This is obtained by matching the rejection region of the UMPBT to that of the fixed design test. For more details please refer to the supplemental information.
}
\usage{
find.umpbt.norm(side = "right", type1 = 0.005, n.obs, null = 0, sigma0 = 1)
}

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
  \item{sigma0}{
a positive numeric; the known standard deviation (sd) in the Z-test

\strong{Default} is \code{1}
}
}

\value{
Returns a real numeric which is the UMPBT point alternative in a Z-test.
}

\references{
MSPRT: main article and supplemental information

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741

Johnson, Valen E., Revised standards for statistical evidence., Proceedings of the National Academy of Sciences, 16, 1945.(Specially it's supplemental file)
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
find.umpbt.norm(n.obs= 60)
}
