\name{design.MSPRT}
\alias{design.MSPRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Designing a modified Sequential Probability Ratio Test (MSPRT)
}
\description{
Given desired values of Type 1 & Type 2 error probabilities and the maximum number of samples (N), this function designs the MSPRT (by finding the 'Termination Threshold' \eqn{\gamma}). \eqn{\gamma} is chosen as the smallest possible value so that the Type 1 error of the MSPRT is maintained at the desired level. This function designs the MSPRT for 3 one sample tests: (1) Test for a binomial proportion, (2) Z-test, and (3) T-test.

By default, this provides the operating characteristics (OC) for the obtained MSPRT at the null hypothesized value. It can also find the same at a user desired point alternative which we can specify through \code{alt.comp}. In general, \code{OC.MSPRT()} can be used to find the OC of a MSPRT at any desired parameter value.

}
\usage{
design.MSPRT(test.type, side, batch.seq, type1 = 0.005, type2 = 0.2,
             null, sigma0 = 1, N.max, alt.comp, repl,
             verbose = T, core.no)
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
  \item{N.max}{
a positive numeric (integer); maximum number of samples that we can afford.
}
  \item{alt.comp}{
missing, \code{TRUE} or a numeric;

If missing, the OC of the MSPRT are obtained only at the \code{null};

If \code{TRUE}, the OC of the MSPRT are computed at the 'fixed design alternative';

If a numeric, it can be any value under the alternative (consistent with \code{"side"}). Then OC of the MSPRT are obtained at this point.
}
  \item{repl}{
a positve numeric (integer); total no. of replications to be used in Monte Carlo method to calculate the OC for an MSPRT;

\strong{Default:} \code{2e+6} for a proportion test; \code{1e+6} for Z or T tests;

should be at least \code{1e+5}.
}
  \item{verbose}{
logical; \code{TRUE} or \code{FALSE};

If \code{TRUE}, returns messages of the current proceedings; otherwise it doesn't;

\strong{Default} is \code{TRUE}.
}
  \item{core.no}{
a numeric;

Number of cores this function can use for carrying out this computation using the parallel computing.

\strong{Default} is \code{1} if there is at most 2 cores, otherwise \code{(number of cores -1)}
}
}

\details{
At the least, we need to provide \code{test.type} and \code{N.max}. For a brief user guide, please refer to the main article and the supplemental information.
}
\value{

If \code{alt.comp} is missing (default), this computes the OC of the obtained MSPRT under the null hypothesis, and returns a list with the following components:
  \item{type1.est}{a numeric in \code{(0,1)}; the Type 1 error probability of the MSPRT.}
  \item{avg.n0}{a positive numeric; the number of samples required on an average by the MSPRT for coming to a decision when the null hypothesis is true.}
  \item{umpbt.alt}{a numeric or a numeric vector of length 2;
  
  For the proportion test, this is usually of length 2. They specify the two points of the UMPBT alternative.
  
  For the Z-test, this is the UMPBT point alternative.
  
  For T-test, this is not returned.
  
}
  
  \item{psi.umpbt}{a numeric in \code{(0,1)};
                   
                  Returned only in case of the proportion test. This denotes the probability of the first component in \code{umpbt.alt}.}
                  
  \item{rej.threshold}{a numeric; the constant value of Wald's rejection threshold.}
  \item{acc.threshold}{a numeric; the constant value of Wald's acceptance threshold.}
  \item{term.thresh}{a positive numeric; denotes the Termination Threshold \eqn{\gamma}.
  }
  
  
If \code{alt.comp} equals \code{TRUE}, this additionally computes the OC of the obtained MSPRT at the 'fixed design alternative'. In this case, the function returns a list with the following components in addition to the previously mentioned components:

  \item{type2.est}{a numeric in \code{(0,1)}; the Type 2 error probability of the obtained MSPRT at the 'fixed design alternative'.}
  \item{avg.n1}{a positive numeric; the number of samples required on an average by the MSPRT for coming to a decision when the 'fixed design alternative' is actually true.}
  \item{alt}{a numeric; the point alternative where the performance is computed.
  
            In this case, this is the 'fixed design alternative'.}
  \item{alt.type2}{a numeric in \code{(0,1)}; the Type 2 error probability of the fixed design test at \code{alt};
                   
                   In this case this is exactly \code{type2}.}


If \code{alt.comp} is numeric, then this computes the list exactly as was for \code{alt.comp = TRUE} with \code{type2.est}, \code{avg.n1}, \code{alt} and \code{alt.type2} now computed at this specified value.

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

## test for a binomial proportion

# can sequentially observe every sample where max available
# sample size is 30
#design.MSPRT(test.type="prop.test", side = "right", null = 0.2, N.max = 30)


# can sequentially observe every fifth sample where max available
# sample size is 30
#design.MSPRT(test.type="prop.test", side = "right", null = 0.2, N.max = 30,
#             batch.seq = seq(5,30, by = 5))
               

               
## Z-test

# can sequentially observe every sample where max available
# sample size is 30
#design.MSPRT(test.type="z.test", side = "right", null = 3, sigma0 = 1.5,
#             N.max = 30)


# can sequentially observe every fifth sample where max available
# sample size is 30
#design.MSPRT(test.type="z.test", side = "right", null = 3, sigma0 = 1.5,
#             N.max = 30, batch.seq = seq(5,30, by = 5))
               


## T-test

# can sequentially observe every sample where max available
# sample size is 30
#design.MSPRT(test.type="t.test", side = "right", null = 3, N.max = 30)


# can sequentially observe every fifth sample where max available
# sample size is 30
#design.MSPRT(test.type="t.test", side = "right", null = 3, N.max = 30,
#             batch.seq = seq(5,30, by = 5))

}
