\name{find.threshold.oneProp}
\alias{find.threshold.oneProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimizing the UMPBT objective function in fixed design one-sample proportion tests
}
\description{
Given \eqn{\delta}, this function finds the difference \code{(optimum value of the objective function - a constant)} in one-sample proportion tests in a fixed design. Notation is as in the supplemental information.
}
\usage{
find.threshold.oneProp(delta, side = "right", n, p0 = 0.5, opt.interval, root)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
a positive numeric; corresponding to this, the UMPBT alternative (UMPBT(\eqn{\delta})) is obtained.
}
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{n}{
a positive numeric (integer); number of samples on which the fixed design one-sample test is based on.
}
  \item{p0}{
a numeric in \code{(0,1)}; the hypothesized value of proportion under the simple null hypothesis.
            
              \strong{Default:} 0.5.
}
  \item{opt.interval}{
a numeric vector of length 2; contains the lower and upper endpoints of an interval to be optimized over.

Endpoints should lie inside \code{(0,1)}.
  
                      \strong{Default:} \code{c(p0,1)} if \code{side = "right"}, and \code{c(0,p0)} if \code{side = "left"}.
}
  \item{root}{
a numeric; the `constant' in the above description.

Given this, the function retunrs \code{(optimum value of the objective function - root)}.
}
}

\details{
Apart from finding the optimum value of the objective function, the argument \code{root} can be used to solve for a \code{delta}. For example, using \code{root=k} we can find a \code{delta} such that the optimized value of the objective function is \code{k}.
}

\value{
If \code{root=k}, this returns a numeric denoting the difference \code{(optimum value of the objective function - k)}.
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
## returns minimum value of the objective function
find.threshold.oneProp(delta= 25, n= 60, p0= 0.2, root= 0)

## returns (minimum value of the objective function -5)
find.threshold.oneProp(delta= 25, n= 60, p0= 0.2, root= 5)
}
