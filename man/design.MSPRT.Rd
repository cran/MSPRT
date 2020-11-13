\name{design.MSPRT}
\alias{design.MSPRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Designing the MSPRT
}
\description{
Given the maximum available sample size and prespecified Type I & II error probabilities, this function designs/obtains the corresponding MSPRT.
}
\usage{
design.MSPRT(test.type, side = "right", theta0, theta1 = T, 
             Type1.target = 0.005, Type2.target = 0.2, 
             N.max, N1.max, N2.max, 
             sigma = 1, sigma1 = 1, sigma2 = 1, 
             batch.size, batch1.size, batch2.size, 
             nReplicate = 1e+06, verbose = T, seed = 1)
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
Numeric. Hypothesized value of effect size (\eqn{\theta_0}) under \eqn{H_0}. 

\strong{Default:} 0.5 in one-sample proportion tests, and 0 for others.
}
  \item{theta1}{
Logical, numeric or list (two components with names \code{'right'} and \code{'left'}).
  
  \itemize{
    \item If \code{FALSE}, no comparison is done under the alternative hypothesis.
    \item If \code{TRUE} (\strong{Default}), comparison is done at the fixed-design alternative effect size (\eqn{\theta_a}).
    \item If numeric, this can only be in case of one-sided tests (that is, \code{side = "right"} or \code{"left"}). The comparison is done at the specified numeric value of the alternative effect size.
    \item If list, this can only be in case of two-sided tests (that is, \code{side = "both"}). The list has to be of the form \code{list("right" =} \eqn{\theta_1}\code{, "left" =} \eqn{\theta_2}). Then the comparison is done at alternative effect sizes \eqn{\theta_1} and \eqn{\theta_2}.
  }

Note: In case of two-sided tests at a given level of significance, there are two effect sizes under \eqn{H_1} (one on the right of \eqn{H_0} and one on the left) that corresponds to the same Type II error probability (or power). This list provides users with the ability where he/she can replace \eqn{\theta_1} and \eqn{\theta_2} by any effect sizes from each side in the form of a list as mentioned above, and can get the designed MSPRT together with its operating characteristics at those effect sizes.
}
  \item{Type1.target}{
Numeric within [0,1]. Prespecified level of Type I error probability.

\strong{Default:} 0.005.

The MSPRT exactly maintains its Type I error probability at this value.
}
  \item{Type2.target}{
Numeric within [0,1]. Prespecified level of Type 2 error probability.

\strong{Default:} 0.2.

The MSPRT approximately maintains its Type II error probability at this value at the corresponding fixed-design alternative (\eqn{\theta_a}).
}
  \item{N.max}{
Positive integer. Maximum available sample size in one-sample tests.
}
  \item{N1.max}{
Positive integer. Maximum available sample size from Group-1 in two-sample tests.
}
  \item{N2.max}{
Positive integer. Maximum available sample size from Group-2 in two-sample tests.
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
  \item{batch.size}{
Integer vector. A vector denoting the number of observations that are planned to be observed at each sequential step in one-sample tests. 
  
\strong{Default:}
  \itemize{
    \item Proportion and z tests: \code{rep(1, N.max)}.
    \item t tests: \code{c(2, rep(1, N.max-1))}.
  }
      
Default values mean the sequential analysis is performed after observing each observation. This corresponds to a sequential MSPRT. If any batch size is more than 1 (or more than 2 in the 1st step for t test) it corresponds to a group sequential MSPRT.

Note: First batch size for t tests needs to be at least 2. The length of batch.size equals to the maximum number of planned sequential analyses.
}
  \item{batch1.size}{
Integer vector. A vector denoting the number of observations that are planned to be observed from Group-1 at each sequential step in two-sample tests. 

\strong{Default:}
  \itemize{
    \item z tests: \code{rep(1, N1.max)}.
    \item t tests: \code{c(2, rep(1, N1.max-1))}.
  }
    
Default values mean the sequential analysis is performed after observing each observation from Group-1.
}
  \item{batch2.size}{
Integer vector. A vector denoting the number of observations that are planned to be observed from Group-2 at each sequential step in two-sample tests. 

\strong{Default:}
  \itemize{
    \item z tests: \code{rep(1, N2.max)}.
    \item t tests: \code{c(2, rep(1, N2.max-1))}.
  }
    
Default values mean the sequential analysis is performed after observing each observation from Group-2.
}
  \item{nReplicate}{
Positive integer. Total number of replications to be used in Monte Carlo simulation for calculating the termination threshold and the operating characteristics of the MSPRT. 

\strong{Default:} 1,000,000.
}
  \item{verbose}{
Logical. If \code{TRUE} (\strong{default}), returns messages of the current proceedings. Otherwise it doesn't.
}
  \item{seed}{
Integer. Random number generating seed. 

\strong{Default:} 1.
}
}
\value{
List. The list has the following named components in case of one-sided one-sample tests:

  \item{TypeI.attained}{
  Numeric in [0,1]. Type I error probability attained by the designed MSPRT.
  }
  
  \item{Type2.attained}{
  Numeric in [0,1]. Type II error probability attained by the designed MSPRT at the specified alternative effect size \code{theta1}. Returned only if \code{theta1} is \code{TRUE} or numeric.
  }
  
  \item{N}{
  List.
  
  \itemize{
    \item If \code{theta1 = FALSE}, the list has one component named \code{H0}. It stores an integer vector of length \code{nReplicate}. This is the vector of sample size required by the MSPRT for each of \code{nReplicate} Monte Carlo simulations under \eqn{H_0}.
    
    \item If \code{theta1} is \code{TRUE} or numeric, the list has two components named \code{H0} and \code{H1}. Each of these stores an integer vector of length \code{nReplicate}. The stored vector under \code{H0} is the same as in \code{theta1 = FALSE}. The \code{H1} component stored the vector of sample size required by the MSPRT for each of \code{nReplicate} Monte Carlo simulations under the specified alternative effect size.
  }
  }
      
  \item{EN}{
  Numeric vector.
  
  \itemize{
    \item If \code{theta1 = FALSE}, the vector is of length 1. It is the number of samples required on average by the MSPRT under \eqn{H_0}.
    
    \item If \code{theta1} is \code{TRUE} or numeric, the vector is of length 2. They are the number of samples required on average by the MSPRT under \eqn{H_0} (first component) and the specified alternative effect size (second component), respectively.
  }
  }
      
  \item{UMPBT or theta.UMPBT}{
  The UMPBT alternative. \code{UMPBT} in case of one-sample proportion test and \code{theta.UMPBT} in case of all the other tests. Their types are the same as their output from \code{UMPBT.alt} function.
  
  Note: Not returned in t tests as it depends on the data.
  }
  
  \item{theta1}{
  Returned only if \code{theta1} is anything but \code{FALSE}. Stores the effect size under \eqn{H_1} where the operating characteristic of the MSPRT is obtained. Of the same type as the argument \code{theta1}.
  }
  
  \item{Type2.fixed.design}{
  Numeric in [0,1]. Type II error probability attained by the fixed design test with sample size \code{N.max} and Type I error probability \code{Type1.target} at the alternative effect size \code{theta1}.
  }
  
  \item{RejectH0.threshold}{
  Positive numeric. Threshold for rejecting \eqn{H_0} in the MSPRT.
  }
  
  \item{RejectH1.threshold}{
  Positive numeric. Threshold for accepting \eqn{H_1} in the MSPRT.
  }
  
  \item{termination.threshold}{
  Positive numeric. Termination threshold of the MSPRT.
  }
  
  In case of one-sided two-sample tests the above components are returned with following modifications:

  \item{N:}{List.
  
    \itemize{
      \item If \code{theta1 = FALSE} the list has one component named \code{H0}.
      
      \item \code{theta1} is \code{TRUE} or numeric, the list has two components named \code{H0} and \code{H1}.
    
      Each of the named components \code{H0} and \code{H1} contains a list with two components named \code{Group1} and \code{Group2}. Each of these contains the same vector corresponding to Group-1 and Group-2. In each of these, it contains the sample size required by the MSPRT in each of \code{nReplicate} Monte Carlo simulations under the respective effect size for the respective group.
    }
  }
  
  \item{EN}{List.
  
    \itemize{
    \item If \code{theta1 = FALSE} the list has one component named \code{H0}.
    
    \item If \code{theta1} is \code{TRUE} or numeric, the list has two components named \code{H0} and \code{H1}.
    }
    
    Each of the named components \code{H0} or \code{H1} contains a list with two components named \code{Group1} and \code{Group2}. In each of these, it contains the sample size required on average by the MSPRT under the respective effect size for the respective group.
    }
  
In case of two-sided tests the above components are returned with following modifications:

    \item{Type2.attained}{
    Numeric vector of length 2 with both elements in [0,1]. The first and second component is the Type II error probability of the MSPRT at the specified alternative effect sizes \code{theta1$right} and \code{theta1$left}, respectively.
    }
  
    \item{N}{
    This is the same as in one-sided tests if \code{theta1} is \code{FALSE}. If \code{theta1} is \code{TRUE} or a two-component list with names \code{right} and \code{left}, this is a list with three components with names \code{H0}, \code{right} and \code{left} instead of a two-component list with names \code{H0} and \code{H1}. Quantities stored under these components are the same as in one-sided tests except the quantities under \code{right} and \code{left} are the same performance of the designed MSPRT at the specified alternative effect sizes \code{theta1$right} and \code{theta1$left}, respectively.
    }
  
    \item{EN}{
    Numeric vector. The same as in one-sided tests if \code{theta1} is \code{FALSE}. If \code{theta1} is \code{TRUE} or a two-component list with names \code{right} and \code{left}, this is a numeric vector of length 3, where the first, second and third components are the average required sample size under \eqn{H_0}, and at the specified alternative effect sizes \code{theta1$right} and \code{theta1$left}, respectively.
    }
  
Additionally, the output list also contains the provided arguments of design.MSPRT, and 
  
    \item{nAnalyses}{
    Positive integer. This is the maximum number of sequential analyses that is planned. This equals to the \code{length(batch.size)} in one-sample tests, and to the \code{length(batch1.size)} and \code{length(batch2.size)} in two-sample tests.
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
#design.MSPRT(test.type = 'oneProp', side = 'right',
#             N.max = 20)

## left-sided
#design.MSPRT(test.type = 'oneProp', side = 'right',
#             N.max = 20)

## two-sided
#design.MSPRT(test.type = 'oneProp', side = 'both',
#             N.max = 20)


##### one-sample z test #####

## right-sided
#design.MSPRT(test.type = 'oneZ', side = 'right',
#             N.max = 20)

## left-sided
#design.MSPRT(test.type = 'oneZ', side = 'right',
#             N.max = 20)

## two-sided
#design.MSPRT(test.type = 'oneZ', side = 'both',
#              N.max = 20)


##### one-sample t test #####

## right-sided
#design.MSPRT(test.type = 'oneT', side = 'right',
#             N.max = 20)

## left-sided
#design.MSPRT(test.type = 'oneT', side = 'right',
#             N.max = 20)

## two-sided
#design.MSPRT(test.type = 'oneT', side = 'both',
#             N.max = 20)


##### two-sample z test #####

## right-sided
#design.MSPRT(test.type = 'twoZ', side = 'right',
#             N1.max = 20, N2.max = 20)

## left-sided
#design.MSPRT(test.type = 'twoZ', side = 'left',
#             N1.max = 20, N2.max = 20)

## two-sided
#design.MSPRT(test.type = 'twoZ', side = 'both',
#             N1.max = 20, N2.max = 20)


##### two-sample t test #####

## right-sided
#design.MSPRT(test.type = 'twoT', side = 'right',
#             N1.max = 20, N2.max = 20)

## left-sided
#design.MSPRT(test.type = 'twoT', side = 'left',
#             N1.max = 20, N2.max = 20)

## two-sided
#design.MSPRT(test.type = 'twoT', side = 'both',
#             N1.max = 20, N2.max = 20)
}
