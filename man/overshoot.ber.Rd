\name{overshoot.ber}
\alias{overshoot.ber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error summary of Wald's SPRT truncated at the maximum available sample size in a proportion test
}
\description{
MSPRT is designed on the assumption that a researcher can afford at most, say, \code{N} samples. This function calculates Type 1 or Type 2 error summaries committed by the Wald's SPRT when it is simply truncated at \code{N} in a proportion test.

It is worth a mention that a case may remain inconcluded due to the truncation. The required sample size for reaching a decision in those cases is \code{N}.
}
\usage{
overshoot.ber(error.type, batch.seq, null, gen.par, alt.LR, alt.psi,
              up, low, N, R, core.no, return.n = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{error.type}{
a character; specifies which of the 2 types of errors need to be accounted for;
            
           \code{"type1"} for Type 1 error;
           
           \code{"type2"} for Type 2 error;
}
  \item{batch.seq}{
a numeric vector;

An increasing sequence until \code{N}. Denotes the sequence of sample sizes where a user will observe data sequentially.

Last element should equal to \code{N}.
}
  \item{null}{
a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple null hypothesis
}
  \item{gen.par}{
a numeric; the value of binomial proportion from which the data needs to be generated from
}
  \item{alt.LR}{
a numeric vector of lenth 2;

this is a vector of the 2-points of the UMPBT alternative; the sequence of weighted likelihood ratios (\eqn{L_n}) is computed in favour of this alternative;

this is same with the output \code{u} from \code{find.umpbt.ber}.
}
  \item{alt.psi}{
a numeric in \code{(0,1)};

denotes the mixing probability corresponding to the first component of \code{alt.LR} in the UMPBT alternative;

this is same with the output \code{psi} from \code{find.umpbt.ber()}.
}
  \item{up}{
a numeric; value of a constant rejection threshold; should be greater than \code{low}.
}
  \item{low}{
a numeric; value of a constant acceptance threshold; should be smaller than \code{up}.
}
  \item{N}{
a positive numeric (integer); number of samples where truncation of Wald's SPRT is required;

in case of a MSPRT, this is the maximum available sample size.
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

\item{count}{a numeric; number of errors of \code{error.type} out of \code{R} replications those are committed}

\item{inconclusive.vec}{a numeric vector; a vector containing values of \eqn{L_N} which remain inconcluded after the truncation}

\item{n.vec}{a numeric vector; a vector of required number of samples; this is of length \code{R}.}


If \code{return.n = FALSE}, the same list with all the components except \code{n.vec} is returned.
}

\references{
Wald, A., Sequential Tests of Statistical Hypotheses. Ann. of Math. Statist., vol. 16, no. 2, 1945, 117-186.
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{

N.max = 30
#overshoot.ber( error.type= "type1", batch.seq= 1:N.max, null= 0.5, gen.par= 0.5,
#               alt.LR= c(0.5,0.55), alt.psi= 0.4, up= 160, low= 0.2, N= N.max,
#               R= 1e+6, core.no= 2, return.n = T)
               
#overshoot.ber( error.type= "type2", batch.seq= 1:N.max, null= 0.5, gen.par= 0.7,
#               alt.LR= c(0.5,0.55), alt.psi= 0.4, up= 160, low= 0.2, N= N.max,
#               R= 1e+6, core.no= 2, return.n = T)

}
