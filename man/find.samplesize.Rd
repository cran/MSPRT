\name{find.samplesize}
\alias{find.samplesize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size required to achieve a higher significance in a fixed design
}
\description{
This function finds the sample size that is required when we decrease the level of significance of a test while still maintaining a desired power at a point alternative. The sample size can be obtained in case of the test for proportion, Z-test and T-test.
}
\usage{
find.samplesize(test.type, N, lower.signif = 0.05, higher.signif = 0.005,
                null.val, side = "right", pow = 0.8, alt, sigma0 = 1,
                n.seq, plot.it = T)
}

\arguments{
  \item{test.type}{
a character; determines the type of test;

             ``\code{prop.test}'' in case of a test for binomial proportion;
             
             ``\code{z.test}'' in case of a Z-test;
             
             ``\code{t.test}'' in case of a T-test.
}
  \item{N}{
a numeric; sample size used at the lower level of significance in a fixed design
}
  \item{lower.signif}{
a numeric in \code{(0,1)}; denotes the lower significance level

\strong{Default} is \code{0.05}.
}
  \item{higher.signif}{
a numeric in \code{(0,1)}; denotes the higher significance level

\strong{Default} is \code{0.005}.
}
  \item{null.val}{
a numeric; 

The hypothesized value of the parameter under the null hypothesis.

The hypothesized parameters are proportion in case of testing a binomial proportion, and  population mean in case of Z and T tests.

\strong{Default:} 0.5 in case of the test for binomial proportion, and 0 for Z and T-tests.
}
  \item{side}{
a character;

direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"}.

\strong{Default} is \code{"right"}.
}
  \item{pow}{
a numeric in \code{(0,1)}; desired level of power at the point alternative

\strong{Default} is 0.8, means \eqn{80} percent power.
}
  \item{alt}{
missing or a numeric; a value of the point alternative where we want to maintain the power;

\strong{Default:} the fixed design alternative at the \code{lower.signif} using \code{N} samples.
}
  \item{sigma0}{
a positive numeric;

Specifies the known standard deviation (sd) in a Z-test;

Only need to be specified in case of a Z-test;

\strong{Default} is 1.
}
  \item{n.seq}{
missing or a numeric vector; final value of the sample size is searched over this vector.

\strong{Default} is \code{N:(4*N)}.
}
  \item{plot.it}{
a logical; if TRUE, returns a plot. Otherwise not.

\strong{Default} is TRUE.
}
}

\value{
Returns a numeric. In a fixed design, this is the sample size that we require to achieve the \code{higher.signif} while still mainting at least \code{pow} amount of power at the \code{alt}.
}
\references{
Main article and supplemental file of MSPRT
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
find.samplesize(test.type = "prop.test", N = 30, null.val = .2)

# In this case, the fixed design alternative at 0.05 is 0.4263. As it seems,
# we need to increase the sample size to 51 to achieve the higher significance
# of 0.005 while still maintaining at least 80% power at 0.4263.


find.samplesize(test.type = "prop.test", N = 30, null.val = .2, side = "left")

# In this case, the fixed design alternative at 0.05 is 0.0516. For testing
# against the left sided alternative, we need to increase the sample size to
# 66 to achieve the higher significance of 0.005 while still maintaining at
# least 80% power at 0.0516.
}
