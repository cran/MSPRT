\name{find.thresold.ber}
\alias{find.thresold.ber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimizing the UMPBT objective function in case of a proportion test
}
\description{
Given \eqn{\delta}, this function finds the difference \code{(optimum value of the objective function - a constant)} in case of a test for binomial proportion. Notation is similar to the supplemental information.
}
\usage{
find.thresold.ber(delta, side = "right", n.obs, p0 = 0.5, opt.interval, root)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{a positive numeric corresponding to which the UMPBT alternative (UMPBT(\eqn{\delta})) is obtained}
  \item{side}{a character;
              
              direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"}.
              
              \strong{Default} is "right"}
  \item{n.obs}{a positive integer; number of samples used in the fixed design test}
  \item{p0}{a numeric in \code{(0,1)}; hypothesized value of the proportion under the simple null hypothesis
            
              \strong{Default: 0.5}}
  \item{opt.interval}{a numeric vector of length 2; contains the lower and upper endpoints of an interval to be optimized over; endpoints should lie inside \code{(0,1)}
  
                      \strong{Default} is \code{c(p0,1)} if \code{side = "right"}, and \code{c(0,p0)} if \code{side = "left"}.}
  \item{root}{a numeric; the \code{constant} in the above description}
}

\details{
Apart from finding the optimum value of the objective function, the argument \code{root} can be used to solve for a \code{delta}. Using \code{root=k} we can find a \code{delta} such that the optimized value of the objective function is \code{k}.
}

\value{
If \code{root=k}, given the \code{delta} this returns a numeric denoting the difference \code{(optimum value of the objective function - k)}
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
## returns minimum value of the objective function
find.thresold.ber(delta= 25, n.obs= 60, p0= 0.2, root= 0)

## returns (minimum value of the objective function -5)
find.thresold.ber(delta= 25, n.obs= 60, p0= 0.2, root= 5)
}

