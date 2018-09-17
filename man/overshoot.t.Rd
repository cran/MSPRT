\name{overshoot.t}
\alias{overshoot.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error summary of Wald's SPRT truncated at the maximum available sample size in a T-test
}
\description{
MSPRT is designed on the assumption that a researcher can at most use, say \code{N}, samples. This function calculates Type 1 or Type 2 error summaries committed by the Wald's SPRT when it is simply truncated at \code{N} in case of a T-test.

It is worth a mention that a case may remain inconcluded due to the truncation. In those cases the required sample size is \code{N}.
}
\usage{
overshoot.t(side, error.type, batch.seq, type1, null, gen.par, up, low,
            N, R, core.no, return.n = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{
a character;

direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"}.
}
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
  \item{type1}{
a numeric in \code{(0,1)};

The probability at which the user wants to control the Type 1 error of the MSPRT.
}
  \item{null}{
a numeric; the hypothesized value of population mean under the simple null hypothesis
}
  \item{gen.par}{
a numeric; the value of population mean from which normal observations need to be generated from.
}
  \item{up}{
a numeric; value of a constant rejection threshold
}
  \item{low}{
a numeric; value of a constant acceptance threshold
}
  \item{N}{
a positive numeric (integer); number of samples where truncation of Wald's SPRT is required
}
  \item{R}{
a positive numeric (integer); number of replications desired in Monte Carlo method; at least \code{1e+5} is required
}
  \item{core.no}{
a numeric;

number of cores this function can use for carrying out required Monte Carlo computation using parallel computing

}
  \item{return.n}{
logical; \code{TRUE} or \code{FALSE}

if \code{TRUE}, corresponding to the cases where Wald's SPRT either accepts or rejects null before truncation, this function will return a vector of number of required samples; otherwise not
}
}

\value{
If \code{return.n = TRUE}, a list with following components is returned:

\item{count}{a numeric; number of errors of \code{error.type} out of \code{R} replications those are committed}

\item{inconclusive.vec}{a numeric vector; a vector containing values of \eqn{L_N} which remain inconcluded after the truncation}

\item{n.vec}{a numeric vector; a vector of required number of samples}


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
#overshoot.t( side="right", error.type= "type1", batch.seq= 1:N.max, null= 0,
#             gen.par= 0, up= 160, low= 0.2, N= N.max,
#             R= 1e+6, core.no= 2, return.n = T)

#overshoot.t( side="right", error.type= "type2", batch.seq= 1:N.max, null= 0,
#             gen.par= 1.5, up= 160, low= 0.2, N= N.max,
#             R= 1e+6, core.no= 2, return.n = T)

}
