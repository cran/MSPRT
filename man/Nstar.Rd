\name{Nstar}
\alias{Nstar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size required to achieve higher significance
}
\description{
Given the sample size that is available at a lower level of significance, this function calculates the sample size that is required for achieving a higher level of significance so that a desired level of Type II error probability is maintained at a desired effect size.
}
\usage{
Nstar(test.type, N, N1, N2, 
      N.increment = 1, N1.increment = 1, N2.increment = 1, 
      lower.signif = 0.05, higher.signif = 0.005, theta0, 
      side = "right", Type2.target = 0.2, theta, 
      sigma = 1, sigma1 = 1, sigma2 = 1, plot.it = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.type}{
Character. Type of test. Currently, the package only allows 

  \itemize{
  
    \item \code{oneProp} for one-sample proportion tests
    \item \code{oneZ} for one-sample z tests
    \item \code{oneT} for one-sample t tests
    \item \code{twoZ} for two-sample z tests
    \item \code{twoT} for two-sample t tests.
  }
}
  \item{N}{
Positive integer. Sample size available at the lower level of significance in one-sample tests.
}
  \item{N1}{
Positive integer. Sample size available from Group-1 at the lower level of significance in two-sample tests.
}
  \item{N2}{
Positive integer. Sample size available from Group-2 at the lower level of significance in two-sample tests.
}
  \item{N.increment}{
Positive integer. Increment in sample size allowed while searching for the sample size that is required for achieving the higher level of significance.
}
  \item{N1.increment}{
Positive integer. Increment in sample size from Group-1 allowed while searching for the sample size that is required for achieving the higher level of significance.
}
  \item{N2.increment}{
Positive integer. Increment in sample size from Group-2 allowed while searching for the sample size that is required for achieving the higher level of significance.
}
  \item{lower.signif}{
Numeric within [0,1]. Lower level of significance. \strong{Default} 0.05.
}
  \item{higher.signif}{
Numeric within [0,1]. Higher level of significance. \strong{Default:} 0.005.
}
  \item{theta0}{
Numeric. Hypothesized value of effect size (\eqn{\theta_0}) under \eqn{H_0}. \strong{Default:} 0.5 in one-sample proportion tests, and 0 for others.
}
  \item{side}{
Character. Direction of the composite alternative hypothesis. \code{right} for \eqn{H_1 : \theta > \theta_0} (\strong{default}), and \code{left} for \eqn{H_1 : \theta < \theta_0}.
}
  \item{Type2.target}{
Numeric within [0,1]. Prespecified level of Type 2 error probability.

\strong{Default:} 0.2.
}
  \item{theta}{
Numeric. Effect size value where \code{Type2.target} Type II error probability is desired at both levels of significance. \strong{Default:} Fixed-design alternative (\eqn{\theta_a}) at the lower level of significance; that is, the effect size where the fixed design test with \eqn{N} samples and level of significance \code{lower.signif} has the Type II error probability \code{Type2.target}.
}
  \item{sigma}{
Positive numeric. Known standard deviation in one-sample z tests. 

\strong{Default:} 1.
}
  \item{sigma1}{
Positive numeric. Known standard deviation for Group-1 in two-sample z tests. 

\strong{Default:} 1.
}
  \item{sigma2}{
Positive numeric. Known standard deviation for Group-2 in two-sample z tests. 

\strong{Default:} 1.
}
  \item{plot.it}{
Logical. If \code{TRUE} (\strong{default}), returns a plot. Otherwise it doesn't.
}
}
\value{
  
  \itemize{
    \item One-sample tests: Numeric. The required sample size.
    
    \item Two-sample tests: Numeric vector of length 2. The first and second components store the sample sizes required respectively from Group 1 and 2 for achieving the higher level of significance.
  }
}
\references{
Pramanik S., Johnson V. E. and Bhattacharya A. (2020+). A Modified Sequential Probability Ratio Test. \href{https://arxiv.org/pdf/1811.08478.pdf}{[Arxiv]}
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}
\examples{

##### one-sample proportion test #####

## right-sided
Nstar(test.type = "oneProp", N = 30)

## left-sided
Nstar(test.type = "oneProp", side = "left", N = 30)


##### one-sample z test #####

## right-sided
Nstar(test.type = "oneZ", N = 30)

## left-sided
Nstar(test.type = "oneZ", side = "left", N = 30)


##### one-sample t test #####

## right-sided
Nstar(test.type = "oneT", N = 30)

## left-sided
Nstar(test.type = "oneT", side = "left", N = 30)


##### two-sample z test #####

## right-sided
Nstar(test.type = "twoZ", N1 = 30, N2 = 30)

## left-sided
Nstar(test.type = "twoZ", side = "left", N1 = 30, N2 = 30)


##### two-sample t test #####

## right-sided
Nstar(test.type = "twoT", N1 = 30, N2 = 30)

## left-sided
Nstar(test.type = "twoT", side = "left", N1 = 30, N2 = 30)
}
