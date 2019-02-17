\name{implement.MSPRT}
\alias{implement.MSPRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implementing a MSPRT
}
\description{
Once we have designed the MSPRT (that is, obtained the termination threshold, \eqn{\gamma}) using \code{design.MSPRT()}, this function implements the MSPRT algorithm for a given data. This is done by sequentially calculating the likelihood ratio(s) or bayes factor(s) (\eqn{L_n}), and then comparing them with the acceptance and rejection thresholds.

This function implements the MSPRT in one-sample proportion tests, and one & two-sample Z & T-tests.
}
\usage{
implement.MSPRT(test.type, obs, obs1, obs2, side, 
                batch.seq, batch1.seq, batch2.seq,
                type1 = 0.005, type2 = 0.2, null, sigma0, term.thresh, 
                N.max, N1.max, N2.max, plot.it = T, verbose = T)
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
  \item{obs}{
a numeric vector; \strong{required only in one-sample tests.}

Denotes the sequentially observed data based on which we want to carry out the null hypothesis significance testing (NHST) using a one-sample test.
}
  \item{obs1}{
a numeric vector; \strong{required only in two-sample tests.}

Denotes the sequentially observed data from Group-1 based on which we want to carry out the NHST using a two-sample test.
}
  \item{obs2}{
a numeric vector; \strong{required only in two-sample tests.}

Denotes the sequentially observed data from Group-2 based on which we want to carry out the NHST using a two-sample test.
}
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{batch.seq}{
a numeric vector; \strong{required only in one-sample tests.}

An increasing sequence of values until \code{N.max}. Denotes the sequence of sample sizes where data is observed sequentially.

Last element should equal to \code{N.max}.

\strong{Default:} In case of one-sample binomial proportion test and one-sample Z-test, this is \code{1:N.max}. The same for one-sample T-test is \code{2:N.max}.
}
  \item{batch1.seq}{
a numeric vector; \strong{required only in two-sample tests.}

An increasing sequence of values until \code{N1.max}. Denotes the sequence of sample sizes where data is observed sequentially from Group-1.

Last element should equal to \code{N1.max}.

\strong{Default:} In case of two-sample Z-test, this is \code{1:N1.max}. The same for two-sample T-test is \code{2:N1.max}.
}
  \item{batch2.seq}{
a numeric vector; \strong{required only in two-sample tests.}

An increasing sequence of values until \code{N2.max}. Denotes the sequence of sample sizes where data is observed sequentially from Group-2.

Last element should equal to \code{N2.max}.

\strong{Default:} In case of two-sample Z-test, this is \code{1:N2.max}. The same for two-sample T-test is \code{2:N2.max}.
}
  \item{type1}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 1 error of the MSPRT.

\strong{Default:} \code{0.005}.
}
  \item{type2}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 2 error of the MSPRT.

\strong{Default:} \code{0.2}.
}
  \item{null}{
a numeric; \strong{required only in one-sample tests.} the hypothesized value of parameter under the null hypothesis.

The hypothesized parameters are proportion in one-sample proportion tests, population mean in one & two-sample Z & T-tests, and difference between the population mean of Group-2 and Group-1 in two-sample Z & T-tests.

In two-sample tests, only \code{null=0} is allowed. This is done automatically. This argument is \strong{ignored} in these cases.

\strong{Default:} 0.5 in one-sample binomial proportion test, and 0 in one-sample Z and T-tests.
}
  \item{sigma0}{
a positive numeric; \strong{required only in one & two-sample Z-tests.} Known population standard deviation in one-sample tests and known common population standard deviation in two-sample tests.

\strong{Default:} 1.
}
  \item{term.thresh}{
a positive numeric; denotes the termination threshold of a MSPRT.

This is determined at the designing step of a MSPRT using \code{design.MSPRT()}.
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
  \item{plot.it}{
logical vector; if \code{TRUE} (\strong{Default}), a comparison plot is returned. Otherwise it's not.
}
  \item{verbose}{
a logical; if \code{TRUE}, returns messages of the current proceedings; otherwise it doesn't.

\strong{Default:} \code{TRUE}.
}
}

\details{
Suppose we want to carry out one of the above tests. To do that, we simply need to follow two steps:

\strong{Step 1: Designing the MSPRT:} Determine the 'Termination threshold' using \code{design.MSPRT()}.

\strong{Step 2: Implementing the MSPRT:} Implement the MSPRT algorithm for a sequentially observed data using \code{implement.MSPRT()}. In Step 2, we need to use the termination threshold obtained from Step 1.
}

\value{
Returns a list with the following components:

\item{decision}{a character; the decision reached based on the provided data.

                Either \code{"reject"} (means `Reject H0') or \code{"accept"} (means `Don't reject H0') or \code{"continue"} (means `continue sampling').}
\item{n}{a numeric (integer); number of samples required for reaching the \code{decision}.}
\item{lhood.ratio}{a numeric vector; a vector of likelihood ratios (\eqn{L_n} 's) in favor of the UMPBT alternative.

This is computed at each element of \code{batch.seq} until the available data or the maximum available sample size.}
\item{rej.threshold}{a numeric; Wald's rejection threshold.}
\item{acc.threshold}{a numeric; Wald's acceptance threshold.}
\item{umpbt.alt}{a numeric or numeric vector; denotes the UMPBT alternative(s).

At each stage, likelihood ratios (\eqn{L_n} 's) are computed in favor of this alternative.}
\item{psi.umpbt}{a numeric; denotes the mixing probability for the UMPBT alternative in one-sample proportion tests.

\strong{Returned only in one-sample proportion tests.}
}
}

\references{
MSPRT: Main article and Supplemental file

Johnson, Valen E., Uniformly most powerful Bayesian tests., Ann. of Stat., 41, (4), 2013, pp. 1716-1741

Johnson, Valen E., Revised standards for statistical evidence., Proceedings of the National Academy of Sciences, 16, 1945.

Daniel J. Benjamin, James O. Berger, Magnus Johannesson, et al. Redefine statistical
significance. Nature Human Behaviour, 2017.
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{

# the termination thresholds are obtained from design.MSPRT()

## one-sample tests

### proportion test

x = c(0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0,
      0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0)

implement.MSPRT(obs = x, test.type = "oneProp", null = 0.2,
                term.thresh = 22.63, N.max = 30)


### Z-test

x = c(4.060319, 5.275465, 3.746557, 7.392921, 5.494262,
      3.769297, 5.731144, 6.107487, 5.863672)
      
implement.MSPRT(obs = x, test.type = "oneZ", null = 3,
                sigma0 = 1.5, term.thresh = 27.856, N.max = 30)


### T-test
x = c(1.738717, 5.076539, 1.116762, 3.105214, 5.567161, 2.095638,
      2.291750, 2.046943, 2.571340, 3.207162, 4.841446, 1.797331)

implement.MSPRT(obs = x, test.type = "oneT", null = 3,
                term.thresh = 32.702, N.max = 30)



## two-sample tests

### Z-test

x = c(4.060319, 5.275465, 3.746557, 7.392921, 5.494262,
      3.769297, 5.731144, 6.107487, 5.863672)
      
implement.MSPRT(obs = x, test.type = "twoZ", null = 3,
                sigma0 = 1.5, term.thresh = , N.max = 30)


### T-test
x = c(1.738717, 5.076539, 1.116762, 3.105214, 5.567161, 2.095638,
      2.291750, 2.046943, 2.571340, 3.207162, 4.841446, 1.797331)

implement.MSPRT(obs = x, test.type = "twoT", null = 3,
                term.thresh = , N.max = 30)

}
