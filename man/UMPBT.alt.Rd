\name{UMPBT.alt}
\alias{UMPBT.alt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
UMPBT alternative
}
\description{
Given a sample size and prespecified Type I & II error probabilities, this function obtains the objective alternative in the Uniformly Most Powerful Bayesian Test (UMPBT).
}
\usage{
UMPBT.alt(test.type, side = "right", theta0, 
          N, N1, N2, Type1 = 0.005, 
          sigma = 1, sigma1 = 1, sigma2 = 1, 
          obs, sd.obs, obs1, obs2, pooled.sd)
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
  \item{side}{
Character. Direction of the composite alternative hypothesis. \code{right} for \eqn{H_1 : \theta > \theta_0} (\strong{default}), and \code{left} for \eqn{H_1 : \theta < \theta_0}.
}
  \item{theta0}{
Numeric. Hypothesized value of effect size (\eqn{\theta_0}) under \eqn{H_0}. \strong{Default:} 0.5 in one-sample proportion tests, and 0 for others.
}
  \item{N}{
Positive integer. Sample size in one-sample tests.
}
  \item{N1}{
Positive integer. Sample size from Group-1 in two-sample tests.
}
  \item{N2}{
Positive integer. Sample size from Group-2 in two-sample tests.
}
  \item{Type1}{
Numeric in [0,1]. Prespecified Type I error probability. \strong{Default:} 0.005.
}
  \item{sigma}{
Positive numeric. Known standard deviation in one-sample z tests. \strong{Default:} 1.
}
  \item{sigma1}{
Positive numeric. Known standard deviation for Group-1 in two-sample z tests. \strong{Default:} 1.
}
  \item{sigma2}{
Positive numeric. Known standard deviation for Group-2 in two-sample z tests. \strong{Default:} 1.
}
  \item{obs}{
Numeric vector. The vector of observations based on which the UMPBT alternative in one-sample t test is determined. Either \code{obs} or \code{sd.obs} is required.
}
  \item{sd.obs}{
Positive numeric. The standard deviation (with divisor n-1) of observations based on which the UMPBT alternative in one-sample t test is determined. Either \code{obs} or \code{sd.obs} is required.
}
  \item{obs1}{
Numeric vector. The vector of observations from Group-1 based on which the UMPBT alternative in two-sample t test is determined. Either both \code{obs1} and \code{obs2}, or \code{pooled.sd} is required.
}
  \item{obs2}{
Numeric vector. The vector of observations from Group-2 based on which the UMPBT alternative in two-sample t test is determined. Either both \code{obs1} and \code{obs2}, or \code{pooled.sd} is required.
}
  \item{pooled.sd}{
Positive numeric. The pooled standard deviation of observations from Group-1 and 2 based on which the UMPBT alternative in two-sample t test is determined. Either both \code{obs1} and \code{obs2}, or \code{pooled.sd} is required.
}
}
\value{
List with two named components \code{theta} and \code{mix.prob} in one-sample proportion test. In this case, the UMPBT alternative is a mixture distribution of two points. \code{theta} contains the two points (effect sizes) and \code{mix.prob} contains their respective mixing probabilities.

Numeric in case of all the other tests. It is the UMPBT alternative effect size.
}
\references{
Johnson, V. E. (2013a). Revised standards for statistical evidence.Proceed-ings of the National Academy of Sciences, 110(48):19313-19317. \href{https://www.stat.tamu.edu/~vjohnson/files/PNAS2013.pdf}{[Article]}

Johnson, V. E. (2013b). Uniformly most powerful Bayesian tests.TheAnnals of Statistics, 41(4):1716-1741. \href{https://projecteuclid.org/euclid.aos/1378386237#info}{[Article]}

Pramanik S., Johnson V. E. and Bhattacharya A. (2020+). A Modified Sequential Probability Ratio Test. \href{https://arxiv.org/pdf/1811.08478.pdf}{[Arxiv]}
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}
\examples{

##### one-sample proportion test #####

## right-sided
UMPBT.alt(test.type = "oneProp", N = 30)

## left-sided
UMPBT.alt(side = "left", test.type = "oneProp", N = 30)


##### one-sample z test #####

## right-sided
UMPBT.alt(test.type = "oneZ", N = 30)

## left-sided
UMPBT.alt(side = "left", test.type = "oneZ", N = 30)


##### one-sample t test #####

## observed data

set.seed(1)
x = rnorm(n = 30, mean = 0, sd = 1.5)
	
## right-sided
UMPBT.alt(test.type = "oneT", N = 30, obs = x)

## left-sided
UMPBT.alt(side = "left", test.type = "oneT", N = 30, obs = x)


##### two-sample z test #####

## right-sided
UMPBT.alt(test.type = "twoZ", N1 = 30, N2 = 30)

## left-sided
UMPBT.alt(side = "left", test.type = "twoZ", N1 = 30, N2 = 30)


##### two-sample t test #####

## observed data

set.seed(1)
x1 = rnorm(n = 30, mean = 0, sd = 1.5)
x2 = rnorm(n = 30, mean = 0, sd = 1.5)

## right-sided
UMPBT.alt(test.type = "twoT", N1 = 30, N2 = 30, 
          obs1 = x1, obs2 = x2)

## left-sided
UMPBT.alt(side = "left", test.type = "twoT", N1 = 30, N2 = 30, 
          obs1 = x1, obs2 = x2)
}
