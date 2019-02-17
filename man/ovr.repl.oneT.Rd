\name{ovr.repl.oneT}
\alias{ovr.repl.oneT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A particular replication step in \code{overshoot.oneT()}
}
\description{
In one-sample T-tests this function simulates a data, computes the bayes factors, and compares with the acceptance and rejection thresholds. \code{overshoot.oneT()} carries out a Monte Carlo method by repeating this function for \code{R} number of times.
}
\usage{
ovr.repl.oneT(side, error.type, batch.seq, type1, null, gen.par,
              up, low, N, seed)
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
ovr.repl.oneT( side = "right", error.type= "type1", batch.seq= 2:N.max,
               type1 = 0.005, null= 0, gen.par= 0, up= 160, low= 0.2, 
               N= N.max, seed= 1)

ovr.repl.oneT( side = "right", error.type= "type2", batch.seq= 2:N.max, 
               type1 = 0.005, null= 0, gen.par= 1.5, up= 160, low= 0.2, 
               N= N.max, seed= 1)

}

