\name{effectiveN.oneProp}
\alias{effectiveN.oneProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating effective maximum sample size to be used in designing the MSPRT in one-sample proportion test
}
\description{
Given a maximum sample size that is planned to use, this function obtains the maximum sample size (\eqn{N}) that is suggested to use in designing the MSPRT for one-sample proportion tests.
}
\usage{
effectiveN.oneProp(N, side = "right", Type1 = 0.005, theta0 = 0.5, 
                   plot.it = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Positive integer. Maximum sample that is intended to use.
}
  \item{side}{
Character. Direction of the composite alternative hypothesis. \code{right} for \eqn{H_1 : \theta > \theta_0} (\strong{default}), and \code{left} for \eqn{H_1 : \theta < \theta_0}.
}
  \item{Type1}{
Numeric in [0,1]. Prespecified Type I error probability. \strong{Default:} 0.005.
}
  \item{theta0}{
Numeric. Hypothesized value of effect size (\eqn{\theta_0}) under \eqn{H_0}. \strong{Default:} 0.5.
}
  \item{plot.it}{
Logical. If \code{TRUE} (\strong{default}), returns a plot. Otherwise it doesn't.
}
}
\value{
Positive integer. This is suggested to use in \code{\link{OCandASN.MSPRT}} as the maximum availeble sample size (\eqn{N}) to design the MSPRT for one-sample proportion tests.
}
\references{
Pramanik S., Johnson V. E. and Bhattacharya A. (2020+). A Modified Sequential Probability Ratio Test. \href{https://arxiv.org/pdf/1811.08478.pdf}{[Arxiv]}
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}
\examples{
effectiveN.oneProp(N = 30)
}
