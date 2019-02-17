\name{overshoot.oneT}
\alias{overshoot.oneT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error summary of the truncated Wald's SPRT in one-sample T-tests
}
\description{
MSPRT is designed on the assumption that we can afford at most, say, \code{N} samples. This function calculates a summary of the Type 1 or Type 2 error committed by the Wald's SPRT when it is simply truncated at \code{N} in one-sample T-tests. This is done using parallel computation.

It is worth a mention that a case may remain inconclusive due to the truncation. The required sample size for reaching that decision is \code{N}.
}
\usage{
overshoot.oneT(side, error.type, batch.seq, type1, null, gen.par,
               up, low, N, R, core.no, return.n = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.
}
  \item{error.type}{
a character; specifies which of the 2 types of errors need to be accounted for.
            
           \code{"type1"} for Type 1 error.
           
           \code{"type2"} for Type 2 error.
}
  \item{batch.seq}{
a numeric vector; an increasing sequence of values until \code{N}. Denotes the sequence of sample sizes where a user will observe data sequentially.

First element should be at least \code{2}. Last element should equal to \code{N}.
}
  \item{type1}{
a numeric in \code{(0,1)}; the probability at which we want to control the Type 1 error of the MSPRT.
}
  \item{null}{
a numeric; the hypothesized value of population mean under the simple null hypothesis.
}
  \item{gen.par}{
a numeric; the value of population mean from which normal observations need to be generated from.
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
#overshoot.oneT( side="right", error.type= "type1", batch.seq= 2:N.max, null= 0,
#                type1 = 0.005, gen.par= 0, up= 160, low= 0.2, N= N.max, R= 1e+6, core.no= 2,
#                return.n = T)

#overshoot.oneT( side="right", error.type= "type2", batch.seq= 2:N.max, null= 0,
#                type1 = 0.005, gen.par= 1.5, up= 160, low= 0.2, N= N.max, R= 1e+6, core.no= 2,
#                return.n = T)

}

