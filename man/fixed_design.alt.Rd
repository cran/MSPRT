\name{fixed_design.alt}
\alias{fixed_design.alt}
\title{
Fixed-design alternative
}
\description{
Given a sample size and prespecified Type I & II error probabilities, this function obtains the fixed-design alternative (\eqn{\theta_a}) for testing the point null hypothesis \eqn{H_0 : \theta = \theta_0}. 
}
\usage{
fixed_design.alt(test.type, side = "right", theta0, 
                 N, N1, N2, Type1 = 0.005, Type2 = 0.2, 
                 sigma = 1, sigma1 = 1, sigma2 = 1)
}
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
  \item{Type2}{
Numeric in [0,1]. Prespecified Type II error probability. \strong{Default:} 0.2.
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
}
\value{
Numeric. The fixed-design alternative effect size (\eqn{\theta_a}).
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
fixed_design.alt(test.type = "oneProp", N = 30)

## left-sided
fixed_design.alt(side = "left", test.type = "oneProp", N = 30)


##### one-sample z test #####

## right-sided
fixed_design.alt(test.type = "oneZ", N = 30)

## left-sided
fixed_design.alt(side = "left", test.type = "oneZ", N = 30)


##### one-sample t test #####

## right-sided
fixed_design.alt(test.type = "oneT", N = 30)

## left-sided
fixed_design.alt(side = "left", test.type = "oneT", N = 30)


##### two-sample z test #####

## right-sided
fixed_design.alt(test.type = "twoZ", N1 = 30, N2 = 30)

## left-sided
fixed_design.alt(side = "left", test.type = "twoZ", N1 = 30, N2 = 30)


##### two-sample t test #####

## right-sided
fixed_design.alt(test.type = "twoT", N1 = 30, N2 = 30)

## left-sided
fixed_design.alt(side = "left", test.type = "twoT", N1 = 30, N2 = 30)
}
