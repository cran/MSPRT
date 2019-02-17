\name{ovr.repl.twoZ}
\alias{ovr.repl.twoZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A particular replication step in \code{overshoot.twoZ()}
}
\description{
In two-sample Z-tests this function simulates a data, computes the likelihood ratios, and compares with the acceptance and rejection thresholds. \code{overshoot.twoZ()} carries out a Monte Carlo method by repeating this function for \code{R} number of times.
}
\usage{
ovr.repl.twoZ(error.type, batch1.seq, batch2.seq, alt.LR, gen.par,
              up, low, N1, N2, seed)
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

N1.max = 30
N2.max = 30

## common sd 1

# both groups have mean 0
ovr.repl.twoZ( error.type= "type1", batch1.seq= 1:N1.max,  batch2.seq= 1:N2.max,
               gen.par= c(0,1), alt.LR = 1, up= 160, low= 0.2,
               N1= N1.max, N2= N2.max, seed= 1)

# Group-1 & 2 have mean -1 and 1, respectively
ovr.repl.twoZ( error.type= "type2", batch1.seq= 1:N1.max,  batch2.seq= 1:N2.max,
               gen.par= c(1,1), alt.LR = 1, up= 160, low= 0.2,
               N1= N1.max, N2= N2.max, seed= 1)

}

