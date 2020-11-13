\name{Type2.fixed_design}
\alias{Type2.fixed_design}
\title{
Type II error probability of fixed design tests
}
\description{
Obtains the Type II error probability of fixed-design tests for testing the point null hypothesis \eqn{H_0 : \theta = \theta_0}.
}
\usage{
Type2.fixed_design(theta, test.type, side = "right", theta0, 
                   N, N1, N2, Type1 = 0.005, 
                   sigma = 1, sigma1 = 1, sigma2 = 1)
}
\arguments{
  \item{theta}{
Numeric. Effect size where the Type II error probability is desired.
}
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
Positive numeric. Known standard deviation in one-sample \eqn{z} tests. \strong{Default:} 1.
}
  \item{sigma1}{
Positive numeric. Known standard deviation for Group-1 in two-sample \eqn{z} tests. \strong{Default:} 1.
}
  \item{sigma2}{
Positive numeric. Known standard deviation for Group-2 in two-sample \eqn{z} tests. \strong{Default:} 1.
}
}
\value{
Numeric in [0,1]. The Type II error probability of the fixed-design test at the specified effect size value \code{theta}.
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
Type2.fixed_design(theta = seq(0, 1, length.out = 10), 
                   test.type = "oneProp", N = 30)

## left-sided
Type2.fixed_design(theta = seq(0, 1, length.out = 10), side = "left", 
                   test.type = "oneProp", N = 30)


##### one-sample z test #####

## right-sided
Type2.fixed_design(theta = seq(0, 1, length.out = 10), 
                   test.type = "oneZ", N = 30)

## left-sided
Type2.fixed_design(theta = seq(-1, 0, length.out = 10), side = "left", 
                   test.type = "oneZ", N = 30)


##### one-sample t test #####

## right-sided
Type2.fixed_design(theta = seq(0, 1, length.out = 10), 
                   test.type = "oneT", N = 30)

## left-sided
Type2.fixed_design(theta = seq(-1, 0, length.out = 10), side = "left", 
                   test.type = "oneT", N = 30)


##### two-sample z test #####

## right-sided
Type2.fixed_design(theta = seq(0, 1, length.out = 10), 
                   test.type = "twoZ", N1 = 30, N2 = 30)

## left-sided
Type2.fixed_design(theta = seq(-1, 0, length.out = 10), side = "left",
                   test.type = "twoZ", N1 = 30, N2 = 30)


##### two-sample t test #####

## right-sided
Type2.fixed_design(theta = seq(0, 1, length.out = 10), 
                   test.type = "twoT", N1 = 30, N2 = 30)

## left-sided
Type2.fixed_design(theta = seq(-1, 0, length.out = 10), side = "left", 
                   test.type = "twoT", N1 = 30, N2 = 30)
}
