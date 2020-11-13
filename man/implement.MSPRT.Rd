\name{implement.MSPRT}
\alias{implement.MSPRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implementing the MSPRT
}
\description{
This function implements the MSPRT for a sequentially observed data.
}
\usage{
implement.MSPRT(obs, obs1, obs2, design.MSPRT.object, 
                termination.threshold, test.type, side = "right",
                theta0, Type1.target = 0.005, Type2.target = 0.2,
                N.max, N1.max, N2.max, 
                sigma = 1, sigma1 = 1, sigma2 = 1,
                batch.size, batch1.size, batch2.size,
                verbose = T, plot.it = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
  Numeric vector. The vector of data in the order they are sequentially observed for one-sample tests.
  Note: Its length can't exceed the length of \code{batch.size}.
}
  \item{obs1}{
  Numeric vector. The vector of data in the order they are sequentially observed from Group-1 for two-sample tests.
  Note: Its length can't exceed the length of \code{batch1.size}.
}
  \item{obs2}{
  Numeric vector. The vector of data in the order they are sequentially observed from Group-2 for two-sample tests.
  Note: Its length can't exceed the length of \code{batch2.size}.
}
  \item{design.MSPRT.object}{
  List. The output returned from \code{\link{design.MSPRT}} corresponding to the MSPRT for which the operating characteristics are desired.
}
  \item{termination.threshold}{
Positive numeric. Termination threshold of the designed MSPRT.
}
  \item{test.type}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{side}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{theta0}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{Type1.target}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{Type2.target}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{N.max}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{N1.max}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{N2.max}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{sigma}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{sigma1}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{sigma2}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{batch.size}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{batch1.size}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{batch2.size}{
Same as in \code{\link{design.MSPRT}}. Not required if \code{design.MSPRT.object} is provided.
}
  \item{verbose}{
Logical. If TRUE (\strong{default}), returns messages of the current proceedings. Otherwise it doesn't.
}
  \item{plot.it}{
  0, 1 or 2 (\strong{default}). 
  \itemize{
    \item if \code{plot.it=0}, no plot is returned.
    \item if \code{plot.it=1}, only the ggplot object required to get a comparison plot is returned, but it's not plotted.
    \item if \code{plot.it=2}, a comparison plot and the corresponding ggplot object is returned.
  }
}
}
\details{
If \code{design.MSPRT.object} is provided, one can only additionally provide \code{nReplicate}, \code{nCore}, \code{verbose} and \code{seed} (Easier option). Otherwise, just like in \code{\link{design.MSPRT}}, all the other arguments together with \code{termination.threshold} (obtained from \code{\link{design.MSPRT}}) needs to be provided adequately.
}
\value{
List. The list has the following named components in case of one-sided one-sample tests:

  \item{n}{
  Positive integer. Number of samples required to reach the decision.
  }
  
  \item{decision}{
  Character. The decision reached. The possibilities are \code{'accept'}, \code{'reject'} and \code{'continue'}. They respectively correspond to accepting \eqn{H_0}, rejecting \eqn{H_0} and continue sampling.
  }
  
  \item{RejectH0.threshold}{
  Positive numeric. Threshold for rejecting \eqn{H_0} in the MSPRT.
  }
  
  \item{RejectH1.threshold}{
  Positive numeric. Threshold for accepting \eqn{H_1} in the MSPRT.
  }
  
  \item{LR}{
  Numeric vector. Vector of weighted likelihood ratios (proportion tests) or likelihood ratios (z tests) or Bayes factor (t tests) that are computed at each step of sequential analysis until either a decision is reached or the maximum available number of samples (\code{N.max} in one-sample tests, or \code{N1.max} and \code{N2.max} in two-sample tests) has been used.
  }
  
  \item{UMPBT alternative}{
  This stores the UMPBT alternative(s) as
  \itemize{
    \item \code{UMPBT} for proportion tests. Of the same type as it is returned by \code{UMPBT.alt} in these tests.
    \item \code{theta.UMPBT} for z and t tests. This is a numeric in case of z tests and a numeric vector in case of t tests. For t tests the UMPBT alternative depends on the data. So the numeric vector returned in this case contains the UMPBT alternative computed at step of sequential analysis and is based on all data observed until that step.
    }
    }
    
  In case of two-sample tests, the \code{n} output above is replaced by \code{n1} and \code{n2}. They are positive integers and refer to the number of samples from Group-1 and 2 required to reach the decision.
  
  In case of two-sided tests at level of significance \eqn{\alpha}, the MSPRT carries out a right and a left sided test simultaneously at level of significance \eqn{\alpha/2}. In this case the outputs are same as above with following changes in components in the returned list:
  
    \item{LR}{
    List. It has two components named \code{right} and \code{left} corresponding to the right and left sided tests of size \eqn{\alpha/2}. Each of these components stores the vector of weighted likelihood ratios (proportion tests) or likelihood ratios (z tests) or Bayes factor (t tests) that are computed at each step of sequential analysis until either a decision is reached or the maximum available number of samples (\code{N.max} in one-sample tests, or \code{N1.max} and \code{N2.max} in two-sample tests) has been used for that sided test.
    }
    
    \item{UMPBT or theta.UMPBT}{List with two components named \code{right} and \code{left} corresponding to the right and left sided tests of size \eqn{\alpha/2}. Each of these contains the UMPBT alternative (of the same type as the output from \code{UMPBT.alt} for the test with respective sides.
    }
}
\references{
Pramanik S., Johnson V. E. and Bhattacharya A. (2020+). A Modified Sequential Probability Ratio Test. \href{https://arxiv.org/pdf/1811.08478.pdf}{[Arxiv]}
}
\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}
\examples{

#################### one-sample proportion test ####################

#### right sided ####
### design
#design.oneprop.right = design.MSPRT(test.type = 'oneProp', side = 'right',
#                                    N.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0.5  # change effect size to experiment
#y = rbinom(20, 1, theta.gen)
#implement.oneprop.right = implement.MSPRT(obs = y, 
#                                          design.MSPRT.object = design.oneprop.right)

#### left sided ####
### design
#design.oneprop.left = design.MSPRT(test.type = 'oneProp', side = 'left',
#                                   N.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0.5  # change effect size to experiment
#y = rbinom(20, 1, theta.gen)
#implement.oneprop.left = implement.MSPRT(obs = y, 
#                                          design.MSPRT.object = design.oneprop.left)

#### both sided ####
### design
#design.oneprop.both = design.MSPRT(test.type = 'oneProp', side = 'both',
#                                   N.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0.5  # change effect size to experiment
#y = rbinom(20, 1, theta.gen)
#implement.oneprop.both = implement.MSPRT(obs = y, 
#                                         design.MSPRT.object = design.oneprop.both)


#################### one-sample z test ####################

#### right sided ####
### design
#design.onez.right = design.MSPRT(test.type = 'oneZ', side = 'right',
#                                 N.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y = rnorm(20, theta.gen, design.onez.right$sigma)
#implement.onez.right = implement.MSPRT(obs = y, 
#                                       design.MSPRT.object = design.onez.right)

#### left sided ####
### design
#design.onez.left = design.MSPRT(test.type = 'oneZ', side = 'left',
#                                N.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y = rnorm(20, theta.gen, design.onez.left$sigma)
#implement.onez.left = implement.MSPRT(obs = y, 
#                                      design.MSPRT.object = design.onez.left)

#### both sided ####
### design
#design.onez.both = design.MSPRT(test.type = 'oneZ', side = 'both',
#                                N.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y = rnorm(20, theta.gen, design.onez.both$sigma)
#implement.onez.both = implement.MSPRT(obs = y, 
#                                      design.MSPRT.object = design.onez.both)


#################### one-sample t test ####################

#### right sided ####
### design
#design.onet.right = design.MSPRT(test.type = 'oneT', side = 'right',
#                                 N.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y = rnorm(20, theta.gen, 1)
#implement.onet.right = implement.MSPRT(obs = y, 
#                                       design.MSPRT.object = design.onet.right)

#### left sided ####
### design
#design.onet.left = design.MSPRT(test.type = 'oneT', side = 'left',
#                                N.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y = rnorm(20, theta.gen, 1)
#implement.onet.left = implement.MSPRT(obs = y, 
#                                      design.MSPRT.object = design.onet.left)

#### both sided ####
### design
#design.onet.both = design.MSPRT(test.type = 'oneT', side = 'both',
#                                N.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y = rnorm(20, theta.gen, 1)
#implement.onet.both = implement.MSPRT(obs = y, 
#                                      design.MSPRT.object = design.onet.both)


#################### two-sample z test ####################

#### right sided ####
### design
#design.twoz.right = design.MSPRT(test.type = 'twoZ', side = 'right',
#                                 N1.max = 20, N2.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y1 = rnorm(20, theta.gen/2, design.twoz.right$sigma1)
#y2 = rnorm(20, -theta.gen/2, design.twoz.right$sigma2)
#implement.twoz.right = implement.MSPRT(obs1 = y1, obs2 = y2,
#                                       design.MSPRT.object = design.twoz.right)

#### left sided ####
### design
#design.twoz.left = design.MSPRT(test.type = 'twoZ', side = 'left',
#                                N1.max = 20, N2.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y1 = rnorm(20, theta.gen/2, design.twoz.left$sigma1)
#y2 = rnorm(20, -theta.gen/2, design.twoz.left$sigma2)
#implement.twoz.left = implement.MSPRT(obs1 = y1, obs2 = y2,
#                                      design.MSPRT.object = design.twoz.left)

#### both sided ####
### design
#design.twoz.both = design.MSPRT(test.type = 'twoZ', side = 'both',
#                                N1.max = 20, N2.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y1 = rnorm(20, theta.gen/2, design.twoz.both$sigma1)
#y2 = rnorm(20, -theta.gen/2, design.twoz.both$sigma2)
#implement.twoz.both = implement.MSPRT(obs1 = y1, obs2 = y2,
#                                      design.MSPRT.object = design.twoz.both)


#################### two-sample t test ####################

#### right sided ####
### design
#design.twot.right = design.MSPRT(test.type = 'twoT', side = 'right',
#                                 N1.max = 20, N2.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y1 = rnorm(20, theta.gen/2, 1)
#y2 = rnorm(20, -theta.gen/2, 1)
#implement.twot.right = implement.MSPRT(obs1 = y1, obs2 = y2,
#                                       design.MSPRT.object = design.twot.right)

#### left sided ####
### design
#design.twot.left = design.MSPRT(test.type = 'twoT', side = 'left',
#                                N1.max = 20, N2.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y1 = rnorm(20, theta.gen/2, 1)
#y2 = rnorm(20, -theta.gen/2, 1)
#implement.twot.left = implement.MSPRT(obs1 = y1, obs2 = y2,
#                                      design.MSPRT.object = design.twot.left)

#### both sided ####
### design
#design.twot.both = design.MSPRT(test.type = 'twoT', side = 'both',
#                                N1.max = 20, N2.max = 20)

### implementation
#set.seed(1)
#theta.gen = 0  # change effect size to experiment
#y1 = rnorm(20, theta.gen/2, 1)
#y2 = rnorm(20, -theta.gen/2, 1)
#implement.twot.both = implement.MSPRT(obs1 = y1, obs2 = y2,
#                                      design.MSPRT.object = design.twot.both)
}
