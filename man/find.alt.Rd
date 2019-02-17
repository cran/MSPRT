\name{find.alt}
\alias{find.alt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding the `fixed design alternative'
}
\description{
This function obtains the point alternative in one-sample binomial proportion tests, one & two-sample Z-tests, and one & two-sample T-tests in a fixed design.
}
\usage{
find.alt(test.type, side = "right", null, n, n1, n2,
         type1 = 0.005, type2 = 0.2, sigma0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.type}{
a character; denotes the type of test;

             ``\code{oneProp}'' for a one-sample binomial proportion test.
             
             ``\code{oneZ}'' for a one-sample Z-test.
             
             ``\code{oneT}'' for a one-sample T-test.
             
             ``\code{twoZ}'' for a two-sample Z-test.
             
             ``\code{twoT}'' for a two-sample T-test.
}
  \item{side}{
a character; direction of the alternative hypothesis H1.

Has to be one of \code{"right"} or \code{"left"}.

\strong{Default:} \code{"right"}.
}
  \item{null}{
a numeric; \strong{required only in one-sample proportion tests and one-sample Z-tests}; hypothesized value of the parameter under the null hypothesis.

The hypothesized parameters are proportion in one-sample proportion tests, population mean in one-sample Z & T-tests, and difference between the population mean of Group-2 and Group-1 in two-sample Z & T-tests.

In one-sample T-tests and two-sample Z & T-tests, only \code{null=0} is allowed. This is done automatically. This argument is \strong{ignored} in these cases.

\strong{Default:} 0.5 in one-sample proportion tests, and 0 in one-sample Z-tests.
}
  \item{n}{
a positive numeric (integer); \strong{required only in one-sample tests}; number of samples on which the fixed design one-sample test is based on.
}
  \item{n1}{
a positive numeric (integer); \strong{required only in two-sample tests}; number of samples from Group-1 on which the fixed design two-sample test is based on.
}
  \item{n2}{
a positive numeric (integer); \strong{required only in two-sample tests}; number of samples from Group-2 on which the fixed design two-sample test is based on.
}
  \item{type1}{
a numeric in \code{(0,1)}; Type 1 error probability of the test.

\strong{Default:} \code{0.005}.
}
  \item{type2}{
a numeric in \code{(0,1)}; Type 2 error probability.

Corresponding to this, the point alternative is obtained.

\strong{Default:} \code{0.2}.
}
  \item{sigma0}{
a positive numeric; \strong{required only in one & two-sample Z-tests.}

The known population standard deviation in one-sample Z-tests.

The known common population standard deviation in two-sample Z-tests.

\strong{Default:} 1.
}
}

\details{
At the `fixed design alternative', the fixed design test of size \code{type1} has Type 2 error \code{type2}.

In one-sample tests, the fixed design is based on \code{n} samples. In two-sample tests, the fixed design is based on \code{n1} & \code{n2} samples from Group 1 & 2, respectively.
}

\value{
Returns a numeric which is the obtained `fixed design alternative'.
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}


\examples{
## finding the alternative in case of a one-sample proportion test
## which provides 80% power against a right and left sided
## alternative, respectively

# default null = 0.5
find.alt(test.type="oneProp", n= 60, type1= 0.005, type2= 0.2)
         
find.alt(test.type="oneProp", side= "left", n= 60,
         type1= 0.005, type2= 0.2)

# null = 0.2         
find.alt(test.type="oneProp", null= 0.2, n= 60, 
         type1= 0.005, type2= 0.2)
         
find.alt(test.type="oneProp", side= "left", null= 0.2,
         n= 60, type1= 0.005, type2= 0.2)


## finding the alternative in case of a one-sample Z-test
## which provides 80% power against a right and left sided
## alternative, respectively
# default sigma0 = 1

# default null = 0
find.alt(test.type="oneZ", n= 60, type1= 0.005, type2= 0.2)

find.alt(test.type="oneZ", side= "left", n= 60, type1= 0.005,
         type2= 0.2)
         
# null = 3
find.alt(test.type="oneZ", null= 3, n= 60, type1= 0.005, type2= 0.2)

find.alt(test.type="oneZ", side= "left", null= 3, n= 60,
         type1= 0.005, type2= 0.2)
         

## finding the alternative in case of a one-sample T-test
## which provides 80% power against a right and left sided
## alternative, respectively

# default null = 0
find.alt(test.type="oneT", n= 60, type1= 0.005, type2= 0.2)
         
find.alt(test.type="oneT", side= "left", n= 60,
         type1= 0.005, type2= 0.2)
         
# null = 3
find.alt(test.type="oneT", null= 3, n= 60, type1= 0.005,
         type2= 0.2)
         
find.alt(test.type="oneT", side= "left", null= 3,
         n= 60, type1= 0.005, type2= 0.2)


## finding the alternative in case of a two-sample Z-test
## which provides 80% power against a right and left sided
## alternative, respectively
# default sigma0 = 1

find.alt(test.type="twoZ", n1= 60, n2= 50, type1= 0.005, type2= 0.2)

find.alt(test.type="twoZ", side= "left", n1= 60, n2= 50, type1= 0.005,
         type2= 0.2)

         
## finding the alternative in case of a two-sample T-test
## which provides 80% power against a right and left sided
## alternative, respectively

find.alt(test.type="twoT", n1= 60, n2= 50, type1= 0.005, type2= 0.2)

find.alt(test.type="twoT", side= "left", n1= 60, n2= 50, type1= 0.005,
         type2= 0.2)


}
