\name{implement.MSPRT}
\alias{implement.MSPRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implementing a modified Sequential Probability Ratio Test (MSPRT)
}
\description{
Once the MSPRT design (or equivalently the termination threshold, \eqn{\gamma}) is obtained using \code{design.MSPRT()}, this function implements the MSPRT algorithm with that \eqn{\gamma} for a given data. This is done by sequentially calculating the likelihood ratio(s) or bayes factor(s) (\eqn{L_n}), and then comparing with the acceptance and rejection thresholds.
}
\usage{
implement.MSPRT(obs, test.type, side, batch.seq, type1 = 0.005, type2 = 0.2,
                null, sigma0, term.thresh, N.max, plot.it = T, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
a numeric vector; denotes the sequentially observed data based on which we want to carry out the null hypothesis significance testing (NHST)

This is a vector of all the observed data (until the present) in the order they are observed.
}
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
  \item{null}{
a numeric; 

The hypothesized value of the parameter under the null hypothesis.

The hypothesized parameters are proportion in case of testing a binomial proportion, and  population mean in case of Z and T tests.

\strong{Default:} 0.5 in case of the test for binomial proportion, and 0 for Z and T-tests.
}
  \item{sigma0}{
a positive numeric;

the known standard deviation (sd) in a Z-test. Only need to be specified when carrying out a Z-test.

\strong{Default} is 1.
}
  \item{term.thresh}{
a positive numeric; the termination threshold;

this is obtained from \code{design.MSPRT()}.
}
  \item{N.max}{
a positive numeric (integer); maximum number of samples that we can afford.
}
  \item{plot.it}{
logical vector; if \code{TRUE} (\strong{Default}) it returns a comparison plot. Otherwise it doesn't.
}
  \item{verbose}{
logical; \code{TRUE} or \code{FALSE};

If \code{TRUE}, returns messages of the current proceedings; otherwise it doesn't;

\strong{Default} is \code{TRUE}.
}
}

\details{
Suppose we want to carry out one of the above one sample tests. We need to follow two steps to carry out a MSPRT for this.

\strong{Step 1: Designing the MSPRT:} Using \code{design.MSPRT()} we determine the 'Termination threshold'.

\strong{Step 2: Implementing the MSPRT:} Using \code{implement.MSPRT()} we implement the MSPRT algorithm for a sequentially observed data. The termination threshold obtained from Step 1 is used here.
}

\value{
Returns a list with the following components:

\item{decision}{a character; the decision reached based on the provided data \code{obs};

                either \code{"reject"} (means reject H0) or \code{"accept"} (means don't reject H0) or \code{"continue"} (means continue sampling)}
\item{n}{a numeric (integer); number of samples required for reaching the \code{decision}}
\item{lhood.ratio}{a numeric vector; a vector of likelihood ratios (\eqn{L_n} 's) in favor of the UMPBT alternative. This is computed at each element of \code{batch.seq} until the available data or \code{N.max}.}
\item{rej.threshold}{a numeric; Wald's rejection threshold}
\item{acc.threshold}{a numeric; Wald's acceptance threshold}
\item{umpbt.alt}{a numeric or numeric vector; denotes the UMPBT alternative(s);

At each stage, likelihood ratios (\eqn{L_n} 's) are computed in favor of this alternative.}
\item{psi.umpbt}{a numeric; denotes the mixing probability for the UMPBT alternative in a proportion test;

returned only in case of a proportion test.}
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

##test for a binomial proportion

x = c(0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0,
      0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0)

implement.MSPRT(obs = x, test.type = "prop.test", null = 0.2,
                term.thresh = 22.63, N.max = 30)


##z-test

x = c(4.060319, 5.275465, 3.746557, 7.392921, 5.494262,
      3.769297, 5.731144, 6.107487, 5.863672)
      
implement.MSPRT(obs = x, test.type = "z.test", null = 3,
                sigma0 = 1.5, term.thresh = 27.856, N.max = 30)


##t-test
x = c(1.738717, 5.076539, 1.116762, 3.105214, 5.567161, 2.095638,
      2.291750, 2.046943, 2.571340, 3.207162, 4.841446, 1.797331)

implement.MSPRT(obs = x, test.type = "t.test", null = 3,
                term.thresh = 32.702, N.max = 30)

}
