\name{check}
\alias{check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequential checking
}
\description{
This compares a sequence of values with two sequences of thresholds, namely upper and lower thresholds. At any step, it rejects if the value gets larger than the upper threshold, accepts if it becomes smaller than the lower threshold, and remains inconclusive if it stays in between until the end. If it stays inconclusive at the last step, this function compares the last value with the termination threshold. If the value is larger than the termination threshold, it rejects; otherwise, it accepts.

As a particular case, this function sequentially checks the likelihood ratios with the upper, lower and termination thresholds according to the MSPRT algorithm (Algorithm 1 in the supplemental file).

This function can be readily adapted to Group sequential designs by specifying \code{batch.seq} (in one-sample tests), and \code{batch1.seq} & \code{batch2.seq} (in two-sample tests) accordingly.
}
\usage{
check(test.type, statistic, upper, lower, 
      batch.seq, batch1.seq, batch2.seq, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.type}{
a character; denotes the type of test;

             ``\code{oneProp}'' for a one-sample proportion test
             
             ``\code{oneZ}'' for a one-sample Z-test
             
             ``\code{oneT}'' for a one-sample T-test
             
             ``\code{twoZ}'' for a two-sample Z-test
             
             ``\code{twoT}'' for a two-sample T-test
}
  \item{statistic}{
a numeric vector; contains a sequence of values which we want to compare.

In particular for the MSPRT, this is a vector of sequentially calculated likelihood ratios (\eqn{L_n}).

It's length should not exceed the length of \code{batch.seq} in one-sample tests, and length of \code{batch1.seq} (or \code{batch2.seq}) in two-sample tests.
}
  \item{upper}{
a numeric vector; the sequence of upper threshold.

It's length should be equal to the length of \code{batch.seq} in one-sample tests, and length of \code{batch1.seq} (or \code{batch2.seq}) in two-sample tests.
}
  \item{lower}{
a numeric vector; the sequence of lower threshold.

It's length should be equal to the length of \code{batch.seq} in one-sample tests, and length of \code{batch1.seq} (or \code{batch2.seq}) in two-sample tests.
}
  \item{batch.seq}{
a numeric vector; \strong{required only in one-sample tests.}
          
          an increasing sequence of sample size values where data are supposed to be observed sequentially
}
  \item{batch1.seq}{
a numeric vector; \strong{required only in two-sample tests.}
          
          an increasing sequence of sample size values where data from Group-1 are supposed to be observed sequentially
}
  \item{batch2.seq}{
a numeric vector; required only in two-sample tests.
          
          an increasing sequence of sample size values where data from Group-2 are supposed to be observed sequentially
}
  \item{threshold}{
a positive numeric; the termination threshold in a MSPRT.
}
}
\details{
Suppose in a one-sample test, we can affrod at most 100 samples and observe the data after every sample. In this case, the \code{batch.seq} will be \code{1:100}.

In another scenario, suppose we observe the data in groups at every 10th sample. In this case, the \code{batch.seq} will be \code{seq( from=10, to=100, by=10)}. So there are at most 10 batches/steps where we can observe the data and compare.

Once the \eqn{L_n}'s, the upper and lower thresholds, and the termination threshold are provided, this function implements the Algorithm 1 in the supplemental file.
}
\value{
In one-sample tests, this returns a list containing the following elements:

  \item{decision}{
a character; the final decision that is made.

Either \code{"accept"}, \code{"reject"} or \code{"continue"}.
}

  \item{n}{
a numeric  (positive integer); number of samples needed for reaching the \code{decision}
}

\item{exit.stage}{
a numeric (positive integer); the step where the decision was reached at
}

In two-sample tests this returns a similar list as above, except now \code{n} is replaced by \code{n1} & \code{n2} having the following descriptions:

  \item{n1}{
a numeric  (positive integer); number of samples needed from Group-1 for reaching the \code{decision}
}
  \item{n2}{
a numeric  (positive integer); number of samples needed from Group-2 for reaching the \code{decision}
}

}


\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}
