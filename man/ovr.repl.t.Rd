\name{ovr.repl.t}
\alias{ovr.repl.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A particular replication step in \code{overshoot.t()}
}
\description{
In case of a T-test, this function generates a data, computes likelihood ratios and compares with acceptance and rejection thresholds. Based on this, \code{overshoot.t()} carries out a Monte Carlo method by repeating this function for, say, \code{R} number of times.
}
\usage{
ovr.repl.t(side, error.type, batch.seq, type1, null, gen.par,
           up, low, N, seed)
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
  \item{seed}{
a positive integer; used in \code{set.seed()} to recreate the simulated data.
}
}

\value{
Returns a list with following components:

\item{incr.count}{either 0 or 1; 1 if and only if an  error of \code{error.type} is made}

\item{inconclusive}{a numeric; the value of \eqn{L_N} if and only if it remains inconclusive after truncating Wald's SPRT at \code{N}; otherwise a numeric of length 0 is returned}

\item{n}{a numeric; if Wald's SPRT either accepts or rejects null hypothesis at or before \code{N}, number of samples required for making that decision is returned; otherwise \code{N} is returned.}
}

\references{
Wald, A., Sequential Tests of Statistical Hypotheses. Ann. of Math. Statist., vol. 16, no. 2, 1945, 117-186.
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{

N.max = 30
ovr.repl.t( side = "right", error.type= "type1", batch.seq= 2:N.max, type1 = 0.005,
            null= 0, gen.par= 0, up= 160, low= 0.2, N= N.max, seed= 1)

ovr.repl.t( side = "right", error.type= "type2", batch.seq= 2:N.max, type1 = 0.005,
            null= 0, gen.par= 1.5, up= 160, low= 0.2, N= N.max, seed= 1)

}
