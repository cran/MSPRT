\name{ovr.repl.oneProp}
\alias{ovr.repl.oneProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A particular replication step in \code{overshoot.oneProp()}
}
\description{
In one-sample proportion tests this function simulates a data, computes the weighted likelihood ratios, and compares with the acceptance and rejection thresholds. \code{overshoot.oneProp()} carries out a Monte Carlo method by repeating this function for \code{R} number of times.
}
\usage{
ovr.repl.oneProp(error.type, batch.seq, null, gen.par, alt.LR, alt.psi,
                 up, low, N, seed)
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
  \item{seed}{
a positive integer; used in \code{set.seed()} to recreate the simulated data.
}
}

\value{
Returns a list with following components:

\item{incr.count}{either 0 or 1; 1 if and only if an error of \code{error.type} is made.}

\item{inconclusive}{a numeric; the value of \eqn{L_N} if and only if it remains inconclusive after truncating Wald's SPRT at \code{N}; otherwise a numeric of length 0 is returned.}

\item{n}{a numeric; number of samples required for reaching the decision. In an inconclusive case, this value is \code{N}.}
}

\references{
Wald, A., Sequential Tests of Statistical Hypotheses. Ann. of Math. Statist., vol. 16, no. 2, 1945, 117-186.
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{

N.max = 30
ovr.repl.oneProp( error.type= "type1", batch.seq= 1:N.max, null= 0.5, gen.par= 0.5,
                  alt.LR= c(0.5,0.55), alt.psi= 0.4, up= 160, low= 0.2, N= N.max,
                  seed= 1)
               
ovr.repl.oneProp( error.type= "type2", batch.seq= 1:N.max, null= 0.5, gen.par= 0.7,
                  alt.LR= c(0.5,0.55), alt.psi= 0.4, up= 160, low= 0.2, N= N.max, 
                  seed= 1)

}
