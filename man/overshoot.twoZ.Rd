\name{overshoot.twoZ}
\alias{overshoot.twoZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error summary of the truncated Wald's SPRT in two-sample Z-tests
}
\description{
MSPRT is designed on the assumption that we can afford at most, say, \code{N} samples. This function calculates a summary of the Type 1 or Type 2 error committed by the Wald's SPRT when it is simply truncated at \code{N} in two-sample Z-tests. This is done using parallel computation.

It is worth a mention that a case may remain inconclusive due to the truncation. The required sample size for reaching that decision is \code{N}.
}
\usage{
overshoot.twoZ(error.type, batch1.seq, batch2.seq, gen.par, 
               alt.LR, up, low, N1, N2, R, core.no, return.n = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{error.type}{
a character; specifies which of the 2 types of errors need to be accounted for.
            
           \code{"type1"} for Type 1 error.
           
           \code{"type2"} for Type 2 error.
}
  \item{batch1.seq}{
a numeric vector; an increasing sequence of values until \code{N1}. Denotes the sequence of sample sizes where data is observed sequentially from Group-1.

Last element should equal be to \code{N1}.
}
  \item{batch2.seq}{
a numeric vector; an increasing sequence of values until \code{N2}. Denotes the sequence of sample sizes where data is observed sequentially from Group-2.

Last element should equal be to \code{N2}.
}
  \item{gen.par}{
a numeric; observations from Group-1 and 2 are generated from the normal distributions with common standard deviation \code{gen.par[2]}, and means \code{-gen.par[1]} and \code{gen.par[1]}, respectively.
}
  \item{alt.LR}{
a numeric; the simple alternative in favor of which the likelihood ratios are calculated sequentially.

The UMPBT point alternative is used in case of a MSPRT. This is same with the output \code{u} from \code{umpbt.twoZ()}.
}
  \item{up}{
a numeric; value of a constant rejection threshold. Should be greater than \code{low}.
}
  \item{low}{
a numeric; value of a constant acceptance threshold. Should be smaller than \code{up}.
}
  \item{N1}{
a positive numeric (integer); maximum available number of samples from Group-1.
}
  \item{N2}{
a positive numeric (integer); maximum available number of samples from Group-2.
}
  \item{R}{
a positive numeric (integer); number of replications desired in the Monte Carlo study.

At least \code{1e+5} is required.
}
  \item{core.no}{
a numeric; number of cores this function can use for carrying out the Monte Carlo study using the parallel computing.
}
  \item{return.n}{
logical; if \code{TRUE}, this returns a vector of sample sizes required for reaching a decision in each of the \code{R} replications.
}
}


\value{
If \code{return.n = TRUE}, a list with following components is returned:

\item{count}{a numeric; number of errors of \code{error.type} those are committed out of \code{R} replications.}

\item{inconclusive.vec}{a numeric vector; a vector containing the likelihood ratio values at \code{N} which remained inconclusive at the last stage.}

\item{n.vec}{a numeric vector; a vector of the required number of samples in each replications. This vector has length \code{R}.}


If \code{return.n = FALSE}, the same list except \code{n.vec} is returned.
}

\references{
Wald, A., Sequential Tests of Statistical Hypotheses. Ann. of Math. Statist., vol. 16, no. 2, 1945, 117-186.
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{

N1.max = 30
N2.max = 30
#overshoot.twoZ( error.type= "type1", batch1.seq= 1:N1.max,  batch2.seq= 1:N2.max,
#                gen.par= c(0,1), alt.LR = 1, up= 160, low= 0.2,
#                N1= N1.max, N2= N2.max, R= 1e+6, core.no= 2, return.n = T)

#overshoot.twoZ( error.type= "type1", batch1.seq= 1:N1.max,  batch2.seq= 1:N2.max,
#                gen.par= c(1,1), alt.LR = 1, up= 160, low= 0.2,
#                N1= N1.max, N2= N2.max, R= 1e+6, core.no= 2, return.n = T)
}


