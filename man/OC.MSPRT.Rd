\name{OC.MSPRT}
\alias{OC.MSPRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Operating characteristics for a MSPRT
}
\description{
This function evaluates the operating characteristics (OC) for a MSPRT at any specified value of the hypothesized parameter. If this specified value lies in the region of the alternative hypothesis, this computes the Type 2 error probability of the MSPRT; otherwise this computes the Type 1 error probability. In both the cases this also computes the number of samples those are required on an average for terminating the sampling and reaching a decision.
}

\usage{
OC.MSPRT(test.type, side, batch.seq, null, term.thresh, theta, sigma0,
         type1 = 0.005, type2 = 0.2, N.max, verbose = T, repl, core.no)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.type}{
a character; determines the type of test;

             ``\code{prop.test}'' in case of a test for binomial proportion;
             
             ``\code{z.test}'' in case of a Z-test;
             
             ``\code{t.test}'' in case of a T-test.
}
  \item{side}{
a character;

direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"}.

\strong{Default} is \code{"right"}.
}
  \item{batch.seq}{
a numeric vector;

An increasing sequence of values until \code{N.max}. Denotes the sequence of sample sizes where a user will observe data sequentially.

Last element should equal to \code{N.max}.

\strong{Default:} In case of the randomized test for binomial proportion and the Z-test,  this is \code{1:N.max}. The same for T-test is \code{2:N.max}.
}
  \item{null}{
a numeric; 

The hypothesized value of the parameter under the null hypothesis.

The hypothesized parameters are proportion in case of testing a binomial proportion, and  population mean in case of Z and T tests.

\strong{Default:} 0.5 in case of the test for binomial proportion, and 0 for Z and T-tests.
}
  \item{term.thresh}{
a positive numeric; denotes the termination threshold of a MSPRT;

this is determined at the designing step of a MSPRT by using \code{design.MSPRT()}.
}
  \item{theta}{
a numeric; denotes the hypothesized parameter value where we want to evaluate the OC for the MSPRT.
}
  \item{sigma0}{
a positive numeric;

the known standard deviation (sd) in a Z-test. Only need to be specified when carrying out a Z-test.

\strong{Default} is 1.
}
  \item{type1}{
a numeric in \code{(0,1)};

The probability at which we want to control the Type 1 error of the MSPRT.

\strong{Default} is \code{0.005}.
}
  \item{type2}{
a numeric in \code{(0,1)};

The probability at which we want to control the Type 2 error of the MSPRT.

\strong{Default} is \code{0.2}.
}
  \item{N.max}{
a positive numeric (integer); maximum number of samples that we can afford.
}
  \item{verbose}{
logical; \code{TRUE} or \code{FALSE};

If \code{TRUE}, returns messages of the current proceedings; otherwise it doesn't;

\strong{Default} is \code{TRUE}.
}
  \item{repl}{
a positve numeric (integer); total no. of replications to be used in Monte Carlo method to calculate the OC for an MSPRT;

\strong{Default:} \code{2e+6} for a proportion test; \code{1e+6} for Z or T tests;

should be at least \code{1e+5}.
}
  \item{core.no}{
a numeric;

Number of cores this function can use for carrying out this computation using the parallel computing.

\strong{Default} is \code{1} if there is at most 2 cores, otherwise \code{(number of cores -1)}
}
}

\details{
For \code{side="right"}, if \code{theta>null} then the Type 2 error is calculated, where as the Type 1 error is calculated when \code{theta<=null}; and vice versa.

To put it simply, if \code{theta} falls in the region of values under the alternative hypothesis, then Type 2 error and \code{avg.n1} are calculated; otherwise Type 1 error and \code{avg.n0} are calculated.
}

\value{
Returns a list with the following components:

\item{type1.est or type2.est}{a numeric in \code{(0,1)}; the Type 1 or Type 2 error probability of the MSPRT evaluated at \code{theta}, respectively.}

\item{avg.n0 or avg.n1}{a positive numeric; the number of samples required on an average by the MSPRT for reaching a decision when \code{theta} is actually true.}

\item{n0.vec or n1.vec}{a numeric vector; denotes a vector of sample sizes required for reaching a decision in each of the \code{repl} replications in the Monte-Carlo study when \code{theta} is actually true.

This is a vector of length \code{repl} .
}

}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{

# the termination thresholds are obtained from design.MSPRT()

## test for a single proportion

#OC.MSPRT(test.type = "prop.test", null = 0.2,
#         term.thresh = 22.63, theta = 0.3, N.max = 30)



## z-test
## finding OC at theta  = 4

#OC.MSPRT(test.type = "z.test", null = 3, sigma0 = 1.5,
#         term.thresh = 27.856, theta = 4, N.max = 30)


## t-test

#OC.MSPRT(test.type = "t.test", null = 3, term.thresh = 32.702,
#         theta = 4, N.max = 30)

}
