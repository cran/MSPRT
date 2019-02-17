\name{overshoot.oneProp}
\alias{overshoot.oneProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error summary of the truncated Wald's SPRT in one-sample proportion tests
}
\description{
MSPRT is designed on the assumption that we can afford at most, say, \code{N} samples. This function calculates a summary of the Type 1 or Type 2 error committed by the Wald's SPRT when it is simply truncated at \code{N} in one-sample proportion tests. This is done using parallel computation.

It is worth a mention that a case may remain inconclusive due to the truncation. The required sample size for reaching that decision is \code{N}.
}
\usage{
overshoot.oneProp(error.type, batch.seq, null, gen.par, alt.LR, alt.psi,
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
a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple null hypothesis.
}
  \item{gen.par}{
a numeric; the value of proportion from which the data needs to be generated from.
}
  \item{alt.LR}{
a numeric vector of lenth 2; this consists of the 2 UMPBT alternative points. The sequence of weighted likelihood ratios (\eqn{L_n}) is computed in favour of this alternative.

This is same with the output \code{u} from \code{umpbt.oneProp()}.
}
  \item{alt.psi}{
a numeric in \code{(0,1)}; the mixing probability corresponding to the first component of \code{alt.LR} in the UMPBT alternative.

This is same with the output \code{psi} from \code{umpbt.oneProp()}.
}
  \item{up}{
a numeric; value of a constant rejection threshold. Should be greater than \code{low}.
}
  \item{low}{
a numeric; value of a constant acceptance threshold. Should be smaller than \code{up}.
}
  \item{N}{
a positive numeric (integer); number of samples where truncation of Wald's SPRT is required.

In a MSPRT, this is the maximum available sample size.
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

N.max = 30
#overshoot.oneProp( error.type= "type1", batch.seq= 1:N.max, null= 0.5,
#                   gen.par= 0.5, alt.LR= c(0.5,0.55), alt.psi= 0.4,
#                   up= 160, low= 0.2, N= N.max, R= 1e+6,
#                   core.no= 2, return.n = T)
               
#overshoot.oneProp( error.type= "type2", batch.seq= 1:N.max, null= 0.5,
#                   gen.par= 0.7, alt.LR= c(0.5,0.55), alt.psi= 0.4,
#                   up= 160, low= 0.2, N= N.max, R= 1e+6, 
#                   core.no= 2, return.n = T)

}

