\name{check}
\alias{check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequential checking
}
\description{
This compares a sequence of values with two sequences of thresholds. If it stays inconclusive after maximum number of comparisons, then it compares with a threshold which is just like the termination threshold in a MSPRT.
}
\usage{
check(test.type, statistic, upper, lower, batch.seq, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.type}{
a character; denotes the type of test;

             ``\code{prop.test}'' in case of a test for binomial proportion
             
             ``\code{z.test}'' in case of a Z-test
             
             ``\code{t.test}'' in case of a T-test
}
  \item{statistic}{
a numeric vector; contains a sequence of values which we want to compare;

In particular for the MSPRT this is a vector of sequentially obtained likelihood ratios (\eqn{L_n}); it's length should not exceed the length of \code{batch.seq}
}
  \item{upper}{
a numeric vector; the sequence of upper threshold; it's length should equal to the length of \code{batch.seq}
}
  \item{lower}{
a numeric vector; the sequence of lower threshold; it's length should equals to the length of \code{batch.seq}
}
  \item{batch.seq}{
a numeric vector; an increasing sequence of sample size values where data are supposed to be observed sequentially
}
  \item{threshold}{
a positive numeric; similar to the termination threshold in a MSPRT
}
}
\details{
Suppose a researcher can affrod at most 100 samples, and can only observe every 10th sample. Then \code{batch.seq} will be \code{seq( from=10, to=100, by=10)}. So there are at most 10 batches/stages where we can observe the data and make compare.

Once the \eqn{L_n}'s, the upper and lower thresholds, and the termination threshold are calculated, this function implements the Algorithm 1 described in the Supplemental file.

Not to mention, any element of the \code{lower} cannot be greater than the corresponding element of the \code{upper}.
}
\value{
Returns a list containing the following elements:

  \item{decision}{
a character; the final decision that is made; either \code{"accept"}, \code{"reject"} or \code{"continue"}
}

  \item{n}{
a numeric; how many samples were needed for coming to the \code{decision}
}
}


\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

