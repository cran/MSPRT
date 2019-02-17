\name{design.MSPRT}
\alias{design.MSPRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Designing a MSPRT
}
\description{
Given the desired values of Type 1 & Type 2 error probabilities and the maximum available number of samples (N), this function designs the MSPRT (by finding the 'Termination Threshold' \eqn{\gamma}). \eqn{\gamma} is the smallest possible value so that the Type 1 error of the MSPRT is maintained at the desired level. This function designs the MSPRT for one-sample proportion tests, one-sample Z and T-tests, and two-sample Z and T-tests.

By default (that is, if \code{alt.comp} is missing or \code{FALSE}) this provides the operating characteristics (OC) for the obtained MSPRT at the null hypothesized value. Otherwise, this also finds the OC at a user desired point alternative. This point alternative is specified through \code{alt.comp}.

In general, \code{OC.MSPRT()} can also be used to find the OC of the MSPRT at any desired parameter value.
}
\usage{
design.MSPRT(test.type, side, batch.seq, batch1.seq, batch2.seq, 
             type1 = 0.005, type2 = 0.2, null, sigma0 = 1, 
             N.max, N1.max, N2.max, alt.comp, repl, 
             verbose = T, core.no)
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
  \item{type1}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 1 error of the MSPRT.

\strong{Default:} \code{0.005}.
}
  \item{type2}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 2 error of the MSPRT.

\strong{Default:} \code{0.2}.
}
  \item{null}{
a numeric; \strong{required only in one-sample proportion and Z-tests.} Denotes value of the  hypothesized parameter under the null hypothesis.

The hypothesized parameters are proportion in one-sample binomial proportion test, population mean in one-sample Z & T-tests, and difference between the population means of Group-2 and Group-1 in two-sample Z & T-tests.

In one-sample T-tests and two-sample tests, only \code{null=0} is allowed. This is done automatically. This argument is \strong{ignored} in these cases.

\strong{Default:} 0.5 in one-sample binomial proportion test, and 0 in one-sample Z-tests.
}
  \item{sigma0}{
a positive numeric; \strong{required only in one & two-sample Z-tests.} Known population standard deviation in one-sample tests and known common population standard deviation in two-sample tests.

\strong{Default:} 1.
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
  \item{alt.comp}{
missing, \code{FALSE} or \code{TRUE}, or a numeric;

If missing or \code{FALSE}, the OC of the MSPRT are obtained only at the \code{null};

If \code{TRUE}, the OC of the MSPRT are computed at the \code{null} and the 'fixed design alternative';

If a numeric, it can be any value under the alternative (consistent with the \code{"side"}). Then OC of the MSPRT are obtained at the \code{null} and this point.
}
  \item{repl}{
a positve numeric (integer); total number of replications to be used in Monte Carlo method to calculate the OC of the MSPRT.

\strong{Default (Recommended):} \code{2e+6} in a one-sample proportion test; \code{1e+6} otherwise. Should be at least \code{1e+5}.
}
  \item{verbose}{
a logical; if \code{TRUE}, returns messages of the current proceedings; otherwise it doesn't.

\strong{Default:} \code{TRUE}.
}
  \item{core.no}{
a numeric; number of cores this function can use for carrying out a parallel computation.

\strong{Default:} \code{1} if there are at most 2 cores, otherwise \code{(number of cores -1)}.
}
}
\details{
In two-sample tests, the hypothesized parameter is (population mean of Group-2 - population mean of Group-1). So \code{null=0} implies that these two means are equal under the null hypothesis; \code{side="right"} implies that the population mean of Group-2 is larger under this alternative hypothesis; \code{side="left"} implies that the population mean of Group-1 is larger under this alternative hypothesis.

For a user guide, please refer to the supplemental information.
}
\value{
Let us first look at the outputs in one-sample tests.

If \code{alt.comp} is missing or \code{FALSE} (default), this computes the OC of the obtained MSPRT under the null hypothesis, and returns a list with the following components:
  \item{type1.est}{a numeric in \code{(0,1)}; the Type 1 error probability of the MSPRT.}
  \item{avg.n0}{a positive numeric; the number of samples required on an average by the MSPRT for coming to a decision when the null hypothesis is true.}
  \item{umpbt.alt}{a numeric or a numeric vector of length 2;
  
  In one-sample proportion tests, this is usually of length 2. They specify the two points of the UMPBT alternative.
  
  In one & two-sample Z-tests, this is the UMPBT point alternative.
  
  In one & two-sample T-tests, this is not returned.
  
}
  
  \item{psi.umpbt}{a numeric in \code{(0,1)}; \strong{returned only in case of one-sample proportion tests.} This denotes the probability of the first component in \code{umpbt.alt}.}
                  
  \item{rej.threshold}{a numeric; the constant value of Wald's rejection threshold.}
  \item{acc.threshold}{a numeric; the constant value of Wald's acceptance threshold.}
  \item{term.thresh}{a positive numeric; denotes the Termination Threshold (\eqn{\gamma}) of a MSPRT.
  }
  
  
If \code{alt.comp} equals \code{TRUE}, this additionally computes the OC of the obtained MSPRT at the 'fixed design alternative'. In this case, the function returns a list with the following components in addition to the previously mentioned components:

  \item{type2.est}{a numeric in \code{(0,1)}; the Type 2 error probability of the obtained MSPRT at the alternative.}
  \item{avg.n1}{a positive numeric; the number of samples required on an average by the MSPRT for coming to a decision when the alternative is true.}
  \item{alt}{a numeric; the point alternative where the performance is computed.
  
            In this case, this is the 'fixed design alternative'.}
  \item{alt.type2}{a numeric in \code{(0,1)}; the Type 2 error probability of the fixed design test at \code{alt}.
                   
                   In this case this is exactly \code{type2}.}


If \code{alt.comp} is numeric, then this computes the list exactly it was for \code{alt.comp = TRUE} with \code{type2.est}, \code{avg.n1}, \code{alt} and \code{alt.type2} now computed at this user specified value.


In two-sample tests we get similar outputs as in one-sample tests above, except \code{avg.n0} and \code{avg.n1} is replaced by \code{avg.n1_0}, \code{avg.n2_0} and \code{avg.n1_1}, \code{avg.n2_1}, respectively.

\item{avg.n1_0}{a positive numeric; the number of samples from Group-1 required on an average by the MSPRT for coming to a decision when the null hypothesis is true.}

\item{avg.n2_0}{a positive numeric; the number of samples from Group-2 required on an average by the MSPRT for coming to a decision when the null hypothesis is true.}

\item{avg.n1_1}{a positive numeric; the number of samples from Group-1 required on an average by the MSPRT for coming to a decision when the alternative hypothesis is true.}

\item{avg.n2_1}{a positive numeric; the number of samples from Group-2 required on an average by the MSPRT for coming to a decision when the alternative hypothesis is true.}

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

## right-sided one-sample proportion test

# can observe data after each sample sequentially. Max available
# sample size is 30
# design.MSPRT(test.type="oneProp", null = 0.2, N.max = 30)


# can observe data at every fifth sample sequentially. Max available
# sample size is 30
# design.MSPRT(test.type="oneProp", null = 0.2, N.max = 30,
#              batch.seq = seq(5,30, by = 5))
               

               
## right-sided one-sample Z-test

# can observe data after each sample sequentially. Max available
# sample size is 30
# design.MSPRT(test.type="oneZ", null = 3, sigma0 = 1.5,
#              N.max = 30)


# can observe data at every fifth sample sequentially. Max available
# sample size is 30
# design.MSPRT(test.type="oneZ", null = 3, sigma0 = 1.5,
#              N.max = 30, batch.seq = seq(5,30, by = 5))
               


## right-sided one-sample T-test

# can observe data after each sample sequentially. Max available
# sample size is 30
# design.MSPRT(test.type="oneT", N.max = 30)


# can observe data at every fifth sample sequentially. Max available
# sample size is 30
# design.MSPRT(test.type="oneT", N.max = 30, batch.seq = seq(5,30, by = 5))


               
## right-sided two-sample Z-test

# can observe data after each sample for each group sequentially.
# Max available sample size is 30 for both groups
# design.MSPRT(test.type="twoZ", sigma0 = 1.5, N1.max = 30, N2.max = 30)


# can observe data at every fifth sample for each group sequentially.
# Max available sample size is 30 for both groups
# design.MSPRT(test.type="twoZ", sigma0 = 1.5,
#              batch1.seq = seq(5,30, by = 5), batch2.seq = seq(5,30, by = 5),
#              N1.max = 30, N2.max = 30)
               

               
## right-sided two-sample T-test

# can observe data after each sample for each group sequentially.
# Max available sample size is 30 for both groups
# design.MSPRT(test.type="twoT", N1.max = 30, N2.max = 30)


# can observe data at every fifth sample for each group sequentially.
# Max available sample size is 30 for both groups
# design.MSPRT(test.type="twoT", N1.max = 30, N2.max = 30,
#              batch1.seq = seq(5,30, by = 5), batch2.seq = seq(5,30, by = 5))
}
