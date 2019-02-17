\name{overshoot.oneZ}
\alias{overshoot.oneZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error summary of the truncated Wald's SPRT in one-sample Z-tests
}
\description{
MSPRT is designed on the assumption that we can afford at most, say, \code{N} samples. This function calculates a summary of the Type 1 or Type 2 error committed by the Wald's SPRT when it is simply truncated at \code{N} in one-sample Z-tests. This is done using parallel computation.

It is worth a mention that a case may remain inconclusive due to the truncation. The required sample size for reaching that decision is \code{N}.
}
\usage{
overshoot.oneZ(error.type, batch.seq, null, gen.par, alt.LR,
               up, low, N, R, core.no, return.n = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{error.type}{
a character; specifies which of the 2 types of errors need to be accounted for.
            
           \code{"type1"} for Type 1 error.
           
           \code{"type2"} for Type 2 error.
}
  \item{batch.seq}{
a numeric vector; an increasing sequence of values until \code{N}. Denotes the sequence of sample sizes where a user will observe data sequentially.

Last element should equal to \code{N}.
}
  \item{null}{
a numeric; the hypothesized value of population mean under the simple null hypothesis.
}
  \item{gen.par}{
a numeric vector of length 2; the first component is the value of the population mean and second component is the known standard deviation (sd).

Observations are generated from a normal distribution with this mean and sd.
}
  \item{alt.LR}{
a numeric; the simple alternative in favor of which the likelihood ratios (\eqn{L_n}) are calculated sequentially.

The UMPBT point alternative is used in case of a MSPRT. This is same with the output \code{u} from \code{umpbt.oneZ()}.
}
  \item{up}{
a numeric; value of a constant rejection threshold; should be greater than \code{low}.
}
  \item{low}{
a numeric; value of a constant acceptance threshold; should be smaller than \code{up}.
}
  \item{N}{
a positive numeric (integer); number of samples where truncation of Wald's SPRT is required;

in a MSPRT, this is the maximum available sample size.
}
  \item{R}{
a positive numeric (integer); number of replications desired in the Monte Carlo study; at least \code{1e+5} is required
}
  \item{core.no}{
a numeric;

number of cores this function can use for carrying out the Monte Carlo study using the parallel computing.
}
  \item{return.n}{
logical;

if \code{TRUE}, this returns a vector of sample sizes required for reaching a decision in each of the \code{R} replications.
}
}

\value{
If \code{return.n = TRUE}, a list with following components is returned:

\item{count}{a numeric; number of errors of \code{error.type} those are committed out of \code{R} replications}

\item{inconclusive.vec}{a numeric vector; a vector containing the values of \eqn{L_N} which remained inconclusive at the last stage}

\item{n.vec}{a numeric vector; a vector of the required number of samples; this is of length \code{R}.}


If \code{return.n = FALSE}, the same list except \code{n.vec} is returned.
}

\references{
Wald, A., Sequential Tests of Statistical Hypotheses. Ann. of Math. Statist., vol. 16, no. 2, 1945, 117-186.
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{

N.max = 30
#overshoot.oneZ( error.type= "type1", batch.seq= 1:N.max, null= 0,
#                gen.par= c(0,1), alt.LR= 1, up= 160, low= 0.2, N= N.max,
#                R= 1e+6, core.no= 2, return.n = T)

#overshoot.oneZ( error.type= "type2", batch.seq= 1:N.max, null= 0,
#                gen.par= c(1.5,1), alt.LR= 1, up= 160, low= 0.2, N= N.max,
#                R= 1e+6, core.no= 2, return.n = T)
}

