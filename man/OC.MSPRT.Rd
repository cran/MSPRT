\name{OC.MSPRT}
\alias{OC.MSPRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Operating characteristics for a MSPRT
}
\description{
This function evaluates the operating characteristics (OC) for a MSPRT at any specified value of the hypothesized parameter. If the specified value lies in the region of the alternative hypothesis, this computes the Type 2 error probability of the MSPRT; otherwise this computes the Type 1 error probability. In both cases, this also computes the number of samples those are required on an average for making a decision.
}

\usage{
OC.MSPRT(test.type, side, batch.seq, batch1.seq, batch2.seq,
         null, term.thresh, theta, sigma0, type1 = 0.005, type2 = 0.2,
         N.max, N1.max, N2.max, verbose = T, repl, core.no)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.type}{
a character; denotes the type of test.

             ``\code{oneProp}'' for a one-sample binomial proportion test.
             
             ``\code{oneZ}'' for a one-sample Z-test.
             
             ``\code{oneT}'' for a one-sample T-test.
             
             ``\code{twoZ}'' for a two-sample Z-test.
             
             ``\code{twoT}'' for a two-sample T-test.
}
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{batch.seq}{
a numeric vector; \strong{required only in one-sample tests.} An increasing sequence of values until \code{N.max}. Denotes the sequence of sample sizes where data is observed sequentially.

Last element should be equal to \code{N.max}.

\strong{Default:} In case of one-sample binomial proportion test and one-sample Z-test, this is \code{1:N.max}. The same for one-sample T-test is \code{2:N.max}.
}
  \item{batch1.seq}{
a numeric vector; \strong{required only in two-sample tests.} An increasing sequence of values until \code{N1.max}. Denotes the sequence of sample sizes where data is observed sequentially from Group-1.

Last element should equal be to \code{N1.max}.

\strong{Default:} In case of two-sample Z-test, this is \code{1:N1.max}. The same for two-sample T-test is \code{2:N1.max}.
}
  \item{batch2.seq}{
a numeric vector; \strong{required only in two-sample tests.} An increasing sequence of values until \code{N2.max}. Denotes the sequence of sample sizes where data is observed sequentially from Group-2.

Last element should equal be to \code{N2.max}.

\strong{Default:} In case of two-sample Z-test, this is \code{1:N2.max}. The same for two-sample T-test is \code{2:N2.max}.
}
  \item{null}{
a numeric; \strong{required only in one-sample tests.} Denotes value of the  hypothesized parameter under the null hypothesis.

The hypothesized parameters are proportion in one-sample binomial proportion test, population mean in one-sample Z & T-tests, and difference between the population means of Group-2 and Group-1 in two-sample Z & T-tests.

In two-sample tests, only \code{null=0} is allowed. This is done automatically. This argument is \strong{ignored} in these cases.

\strong{Default:} 0.5 in one-sample binomial proportion test, and 0 in one-sample Z-tests.
}
  \item{term.thresh}{
a positive numeric; denotes the termination threshold of a MSPRT.

This is determined at the designing step of a MSPRT using \code{design.MSPRT()}.
}
  \item{theta}{
a numeric; the value of hypothesized parameter where we want to evaluate the OC for the MSPRT.
}
  \item{sigma0}{
a positive numeric; \strong{required only in one & two-sample Z-tests.} Known population standard deviation in one-sample tests and known common population standard deviation in two-sample tests.

\strong{Default:} 1.
}
  \item{type1}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 1 error of the MSPRT.

\strong{Default:} \code{0.005}.
}
  \item{type2}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 2 error of the MSPRT.

\strong{Default:} \code{0.2}.
}
  \item{N.max}{
a positive numeric (integer); \strong{required only in one-sample tests.} Maximum number of samples that we can afford in the one-sample test.
}
  \item{N1.max}{
a positive numeric (integer); \strong{required only in two-sample tests.} Maximum number of samples from Group-1 that we can afford in the two-sample test.
}
  \item{N2.max}{
a positive numeric (integer); \strong{required only in two-sample tests.} Maximum number of samples from Group-2 that we can afford in the two-sample test.
}
  \item{verbose}{
a logical; if \code{TRUE}, returns messages of the current proceedings; otherwise it doesn't.

\strong{Default:} \code{TRUE}.
}
  \item{repl}{
a positve numeric (integer); total number of replications to be used in Monte Carlo method to calculate the OC of the MSPRT.

\strong{Default (Recommended):} \code{2e+6} in a one-sample proportion test; \code{1e+6} otherwise. Should be at least \code{1e+5}.
}
  \item{core.no}{
a numeric; number of cores this function can use for carrying out a parallel computation.

\strong{Default:} \code{1} if there are at most 2 cores, otherwise \code{(number of cores -1)}.
}
}

\details{
For \code{side="right"}, if \code{theta>null} then the Type 2 error is calculated, where as the Type 1 error is calculated when \code{theta<=null}; and vice versa.

To put it simply, if \code{theta} falls under the alternative hypothesis, then Type 2 error and \code{avg.n1} (in one-sample tests) or \code{avg.n1_1} & \code{avg.n2_1} (in two-sample tests) are calculated; otherwise Type 1 error and \code{avg.n0} (in one-sample tests) or \code{avg.n1_0} & \code{avg.n2_0} (in two-sample tests) are calculated.
}

\value{
Returns a list with the following components:

In one-sample tests:

\item{type1.est or type2.est}{a numeric in \code{(0,1)}; the Type 1 or Type 2 error probability of the MSPRT evaluated at \code{theta}, respectively.}

\item{avg.n0 or avg.n1}{a positive numeric; the number of samples required on an average by the MSPRT for reaching a decision when \code{theta} is actually true.}

\item{n0.vec or n1.vec}{a numeric vector; denotes a vector of sample sizes required for reaching a decision in each of the \code{repl} replications in the Monte-Carlo study when \code{theta} is actually true.

This is a vector of length \code{repl} .
}


In two-sample tests:

\item{type1.est or type2.est}{a numeric in \code{(0,1)}; the Type 1 or Type 2 error probability of the MSPRT evaluated at \code{theta}, respectively.}

\item{avg.n1_0 & avg.n2_0, or avg.n1_1 & avg.n2_1}{a positive numeric; the number of samples required on an average by the MSPRT from Group 1 and 2 for reaching a decision when \code{theta} is actually true.}

\item{n1_0.vec & n2_0.vec, or n1_1.vec & n2_1.vec}{a numeric vector; denotes a vector of sample sizes required from Group 1 and 2 for reaching a decision in each of the \code{repl} replications in the Monte-Carlo study when \code{theta} is actually true.

This is a vector of length \code{repl} .
}

}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{

# the termination thresholds are obtained from design.MSPRT()

## One-sample proportion test
## finding OC at theta  = 0.3

# OC.MSPRT(test.type = "oneProp", null = 0.2,
#          term.thresh = 22.63, theta = 0.3, N.max = 30)



## One-sample Z-test
## finding OC at theta  = 4

# OC.MSPRT(test.type = "oneZ", null = 3, sigma0 = 1.5,
#          term.thresh = 27.856, theta = 4, N.max = 30)


## One-sample T-test
## finding OC at theta  = 4

# OC.MSPRT(test.type = "oneT", term.thresh = 32.702,
#          null = 3, theta = 4, N.max = 30)


## Two-sample Z-test
## finding OC at theta  = 1

# OC.MSPRT(test.type = "twoZ", sigma0 = 1.5, term.thresh = 27.163,
#          theta = 1, N1.max = 30, N2.max = 30)


## Two-sample T-test
## finding OC at theta  = 1

# OC.MSPRT(test.type = "twoT", term.thresh = 32.773,
#          theta = 1, N1.max = 30, N2.max = 30)

}

