\name{find.samplesize}
\alias{find.samplesize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size required to achieve a higher significance for one-sample tests in a fixed design
}
\description{
This function finds the sample size that is required to maintain a desired power at a point alternative when we decrease the level of significance of a fixed design one-sample test. This can be calculated in one-sample proportion tests, and one-sample Z & T-tests.
}
\usage{
find.samplesize(test.type, N, lower.signif = 0.05, higher.signif = 0.005,
                null, side = "right", pow = 0.8, alt, sigma0 = 1,
                n.seq, verbose=T, plot.it = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.type}{
a character; denotes the type of test;

             ``\code{oneProp}'' for a one-sample proportion tests.
             
             ``\code{oneZ}'' for a one-sample Z-tests.
             
             ``\code{oneT}'' for a one-sample T-tests.
}
  \item{N}{
a positive numeric (integer); \strong{required only in one-sample tests}; number of samples on which the fixed design one-sample test is based on.
}
  \item{lower.signif}{
a numeric in \code{(0,1)}; denotes the lower level of significance.

\strong{Default:} \code{0.05}.
}
  \item{higher.signif}{
a numeric in \code{(0,1)}; denotes the higher level of significance.

\strong{Default:} \code{0.005}.
}
  \item{null}{
a numeric; denotes the hypothesized value of the parameter under the null hypothesis.

The hypothesized parameters are proportion in one-sample proportion tests, and population mean in one-sample Z & T-tests.

\strong{Default:} 0.5 in one-sample proportion tests, and 0 in one-sample Z and T-tests.
}
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{pow}{
a numeric in \code{(0,1)}; desired level of power at the point alternative \code{alt}.

\strong{Default:} 0.8. This means \eqn{80} percent power.
}
  \item{alt}{
missing or a numeric; value of the point alternative where we want to maintain the power \code{pow}.

\strong{Default:} the `fixed design alternative' at the \code{lower.signif} using \code{N} samples.
}
  \item{sigma0}{
a positive numeric; \strong{required only in one-sample Z-tests.} known population standard deviation.

\strong{Default:} 1.
}
  \item{n.seq}{
missing or a numeric vector; the final value of increased sample size is searched over this vector.

\strong{Default:} \code{N:(4*N)}.
}
  \item{verbose}{
logical; if \code{TRUE}, returns a message of the achieved results; otherwise it doesn't.

\strong{Default:} \code{TRUE}.
}
  \item{plot.it}{
a logical; if \code{TRUE}, returns a plot; otherwise it doesn't.

\strong{Default:} \code{TRUE}.
}
}

\value{
Returns a numeric. In a fixed design for the specified one-sample test, this is the sample size that we require to achieve \code{higher.signif} while still mainting at least \code{pow} amount of power at \code{alt}.
}
\references{
Main article and supplemental file of MSPRT
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
find.samplesize(test.type = "oneProp", N = 30, null = .2)

# In this case, the fixed design alternative at 0.05 is 0.4263. As it seems,
# we need to increase the sample size to 51 to achieve the higher significance
# of 0.005 while still maintaining at least 80% power at 0.4263.


find.samplesize(test.type = "oneProp", N = 30, null = .2, side = "left")

# In this case, the fixed design alternative at 0.05 is 0.0516. For testing
# against the left sided alternative, we need to increase the sample size to
# 66 to achieve the higher significance of 0.005 while still maintaining at
# least 80% power at 0.0516.

find.samplesize(test.type = "oneZ", N = 30)
find.samplesize(test.type = "oneZ", N = 30, side = "left")

find.samplesize(test.type = "oneT", N = 30)
find.samplesize(test.type = "oneT", N = 30, side = "left")
}
