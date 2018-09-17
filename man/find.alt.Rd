\name{find.alt}
\alias{find.alt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding the `fixed design alternative'
}
\description{
This function finds the `fixed design alternative' in case of a test for binomial proportion, Z and T-test.
}
\usage{
find.alt(test.type, side = "right", null, size, type1 = 0.005, type2 = 0.2, sigma0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.type}{
a character;

            ``\code{prop.test}'' in case of a test for binomial proportion
             
             ``\code{z.test}'' in case of a Z-test
             
             ``\code{t.test}'' in case of a T-test
}
  \item{side}{
a character; direction of the alternative hypothesis H1; has to be one of \code{"right"} or \code{"left"}
}
  \item{null}{
a numeric; denotes the hypothesized value under the null hypothesis;

The hypothesized parameters are proportion in case of testing a binomial proportion, and  population mean in case of Z and T tests.

has to be in \code{(0,1)} in case of a proportion test, but can be any real value in case of a Z or T Test
}
  \item{size}{
a positive numeric (integer); number of samples the fixed design test is based on
}
  \item{type1}{
a numeric in \code{(0,1)}; desired Type 1 error probability of the fixed design test
}
  \item{type2}{
a numeric in \code{(0,1)}; desired Type 2 error probability corresponding to which we want to obtain the `fixed design alternative'
}
  \item{sigma0}{
a positive numeric; denotes the known standard deviation in a Z-test;

need to be specified only if \code{test.type = "z.test"}
}
}


\value{
Returns a numeric which is the obtained `fixed design alternative'.
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}


\examples{
## finding the alternative in case of a proportion test which
## provides 80% power against a right and left sided
## alternative, respectively

# default null = 0.5
find.alt(test.type="prop.test", side= "right", size= 60,
         type1= 0.005, type2= 0.2)
         
find.alt(test.type="prop.test", side= "left", size= 60,
         type1= 0.005, type2= 0.2)

# null = 0.2         
find.alt(test.type="prop.test", side= "right", null= 0.2,
         size= 60, type1= 0.005, type2= 0.2)
         
find.alt(test.type="prop.test", side= "left", null= 0.2,
         size= 60, type1= 0.005, type2= 0.2)


## finding the alternative in case of a Z-test which
## provides 80% power against a right and left sided
## alternative, respectively

# default null = 0
find.alt(test.type="z.test", side= "right",
         size= 60, type1= 0.005, type2= 0.2, sigma0= 1)

find.alt(test.type="z.test", side= "left",
         size= 60, type1= 0.005, type2= 0.2, sigma0= 1)
         
# null = 3
find.alt(test.type="z.test", side= "right", null= 3,
         size= 60, type1= 0.005, type2= 0.2, sigma0= 1)

find.alt(test.type="z.test", side= "left", null= 3,
         size= 60, type1= 0.005, type2= 0.2, sigma0= 1)
         

## finding the alternative in case of a T-test which
## provides 80% power against a right and left sided
## alternative, respectively

# default null = 0
find.alt(test.type="t.test", side= "right", size= 60,
         type1= 0.005, type2= 0.2)
         
find.alt(test.type="t.test", side= "left", size= 60,
         type1= 0.005, type2= 0.2)
         
# null = 3
find.alt(test.type="t.test", side= "right", null= 3,
         size= 60, type1= 0.005, type2= 0.2)
         
find.alt(test.type="t.test", side= "left", null= 3,
         size= 60, type1= 0.005, type2= 0.2)

}
