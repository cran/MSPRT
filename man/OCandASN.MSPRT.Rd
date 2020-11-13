\name{OCandASN.MSPRT}
\alias{OCandASN.MSPRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Operating characteristics (OC) and Average Sample Number (ASN) of a designed MSPRT
}
\description{
This function obtains the operating characteristics, that is the probability of accepting \eqn{H_0} and the sample size required on average for reaching a decision, for a designed MSPRT at the specified effect size(s).
}
\usage{
OCandASN.MSPRT(theta, design.MSPRT.object, 
               termination.threshold, test.type, side = "right", 
               theta0, Type1.target = 0.005, Type2.target = 0.2,
               N.max, N1.max, N2.max,
               sigma = 1, sigma1 = 1, sigma2 = 1, 
               batch.size, batch1.size, batch2.size, 
               nReplicate = 1e+06, nCore = max(1, detectCores() - 1),
               verbose = T, seed = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
Numeric vector. Vector of effect size(s) where the operating characteristics of the MSPRT is desired.
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
  \item{nReplicate}{
Positive integer. Total number of replications to be used in Monte Carlo simulation for calculating the termination threshold and the operating characteristics of the MSPRT. 

\strong{Default:} 1,000,000.
}
  \item{verbose}{
Logical. If \code{TRUE} (\strong{default}), returns messages of the current proceedings. Otherwise it doesn't.
}
  \item{nCore}{
Positive integer. Total number of cores available for computation. Can be anything \eqn{\ge 1}. 

  \strong{Default:} \code{detectCores() - 1}. That is, 1 less than the total number of available cores.
}
  \item{seed}{
Integer. Random number generating seed. 

\strong{Default:} 1.
}
}
\details{
If \code{design.MSPRT.object} is provided, one can only additionally provide \code{nReplicate}, \code{nCore}, \code{verbose} and \code{seed} (Easier option). Otherwise, just like in \code{\link{design.MSPRT}}, all the other arguments together with \code{termination.threshold} (obtained from \code{\link{design.MSPRT}}) needs to be provided adequately.
}
\value{
  Data frame.
  
  \itemize{
    \item One-sample tests: The data frame has 3 columns named \code{theta}, \code{acceptH0.prob} and \code{EN}, and the number of rows equals to the number of effect sizes (length of \code{theta}) where the operating characteristics are evaluated. Each row corresponds to a particular value of theta (effect size). The columns respectively contain the value of a particular theta (effect size), and the probability of accepting the $H_0$ and the average sample size required by the MSPRT for reaching a decision thereat.
    
    \item Two-sample tests: The data frame has 4 columns named \code{theta}, \code{acceptH0.prob}, \code{EN1} and \code{EN2}, and the number of rows equals to the number of effect sizes (length of \code{theta}) where the operating characteristics are evaluated. Each row corresponds to a particular value of theta (effect size). The columns respectively contain the value of a particular theta (effect size), and the probability of accepting the \eqn{H_0} at that effect size, and the average sample size from Group-1 & 2 that is required by the MSPRT for reaching a decision thereat.
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
#design.oneprop.right <- design.MSPRT(test.type = 'oneProp', side = 'right',
#                                     N.max = 20)

### OC and ASN
#OC.oneprop.right <- OCandASN.MSPRT(theta = seq(design.oneprop.right$theta0, 1,
#                                               length.out = 3),
#                                   design.MSPRT.object = design.oneprop.right)

#### left sided ####
### design
#design.oneprop.left = design.MSPRT(test.type = 'oneProp', side = 'left',
#                                   N.max = 20)

### OC and ASN
#OC.oneprop.left = OCandASN.MSPRT(theta = seq(0, design.oneprop.right$theta0,
#                                             length.out = 3),
#                                 design.MSPRT.object = design.oneprop.left)

#### both sided ####
### design
#design.oneprop.both = design.MSPRT(test.type = 'oneProp', side = 'both',
#                                   N.max = 20)

### OC and ASN
#OC.oneprop.both = OCandASN.MSPRT(theta = seq(0, 1, length.out = 3),
#                                 design.MSPRT.object = design.oneprop.both)


#################### one-sample z test ####################

#### right sided ####
### design
#design.onez.right = design.MSPRT(test.type = 'oneZ', side = 'right',
#                                 N.max = 20)

### OC and ASN
#OC.onez.right = OCandASN.MSPRT(theta = seq(design.onez.right$theta0,
#                                           design.onez.right$theta0 + 3*design.onez.right$sigma,
#                                           length.out = 3),
#                               design.MSPRT.object = design.onez.right)

#### left sided ####
### design
#design.onez.left = design.MSPRT(test.type = 'oneZ', side = 'left',
#                                N.max = 20)

### OC and ASN
#OC.onez.left = OCandASN.MSPRT(theta = seq(design.onez.left$theta0 - 3*design.onez.left$sigma,
#                                          design.onez.left$theta0,
#                                          length.out = 3),
#                              design.MSPRT.object = design.onez.left)

#### both sided ####
### design
#design.onez.both = design.MSPRT(test.type = 'oneZ', side = 'both',
#                                N.max = 20)

### OC and ASN
#OC.onez.both = OCandASN.MSPRT(theta = seq(design.onez.both$theta0 - 3*design.onez.both$sigma,
#                                          design.onez.both$theta0 + 3*design.onez.both$sigma,
#                                          length.out = 3),
#                              design.MSPRT.object = design.onez.both)


#################### one-sample t test ####################

#### right sided ####
### design
#design.onet.right = design.MSPRT(test.type = 'oneT', side = 'right',
#                                 N.max = 20)

### OC and ASN
#OC.onet.right = OCandASN.MSPRT(theta = seq(design.onet.right$theta0, 1,
#                                           length.out = 3),
#                               design.MSPRT.object = design.onet.right)

#### left sided ####
### design
#design.onet.left = design.MSPRT(test.type = 'oneT', side = 'left',
#                                N.max = 20)

### OC and ASN
#OC.onet.left = OCandASN.MSPRT(theta = seq(-1, design.onet.left$theta0,
#                                          length.out = 3),
#                              design.MSPRT.object = design.onet.left)

#### both sided ####
### design
#design.onet.both = design.MSPRT(test.type = 'oneT', side = 'both',
#                                N.max = 20)

### OC and ASN
#OC.onet.both = OCandASN.MSPRT(theta = seq(-1, 1, length.out = 3),
#                              design.MSPRT.object = design.onet.both)


#################### two-sample z test ####################

#### right sided ####
### design
#design.twoz.right = design.MSPRT(test.type = 'twoZ', side = 'right',
#                                 N1.max = 20, N2.max = 20)

### OC and ASN
#OC.twoz.right = OCandASN.MSPRT(theta = seq(design.twoz.right$theta0, 
#                                           design.twoz.right$theta0 + 2,
#                                           length.out = 3),
#                               design.MSPRT.object = design.twoz.right)

#### left sided ####
### design
#design.twoz.left = design.MSPRT(test.type = 'twoZ', side = 'left',
#                                N1.max = 20, N2.max = 20)

### OC and ASN
#OC.twoz.left = OCandASN.MSPRT(theta = seq(design.twoz.left$theta0 - 2,
#                                          design.twoz.left$theta0,
#                                          length.out = 3),
#                              design.MSPRT.object = design.twoz.left)

#### both sided ####
### design
#design.twoz.both = design.MSPRT(test.type = 'twoZ', side = 'both',
#                                N1.max = 20, N2.max = 20)

### OC and ASN
#OC.twoz.both = OCandASN.MSPRT(theta = seq(design.twoz.both$theta0 - 2,
#                                          design.twoz.both$theta0 + 2,
#                                          length.out = 3),
#                             design.MSPRT.object = design.twoz.both)


#################### two-sample t test ####################

#### right sided ####
### design
#design.twot.right = design.MSPRT(test.type = 'twoT', side = 'right',
#                                 N1.max = 20, N2.max = 20)

### OC and ASN
#OC.twot.right = OCandASN.MSPRT(theta = seq(design.twot.right$theta0, 
#                                           design.twot.right$theta0 + 2,
#                                           length.out = 3),
#                               design.MSPRT.object = design.twot.right)

#### left sided ####
### design
#design.twot.left = design.MSPRT(test.type = 'twoT', side = 'left',
#                                N1.max = 20, N2.max = 20)

### OC and ASN
#OC.twot.left = OCandASN.MSPRT(theta = seq(design.twot.left$theta0 - 2,
#                                          design.twot.left$theta0,
#                                          length.out = 3),
#                              design.MSPRT.object = design.twot.left)

#### both sided ####
### design
#design.twot.both = design.MSPRT(test.type = 'twoT', side = 'both',
#                                N1.max = 20, N2.max = 20)

### OC and ASN
#OC.twot.both = OCandASN.MSPRT(theta = seq(design.twot.both$theta0 - 2,
#                                          design.twot.both$theta0 + 2,
#                                          length.out = 3),
#                              design.MSPRT.object = design.twot.both)
}
