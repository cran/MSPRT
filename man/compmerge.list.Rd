\name{compmerge.list}
\alias{compmerge.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merging two lists componentwise
}
\description{
Suppose two lists have same number of components. This function creates a new list by merging (or concatenating) them componentwise.
}
\usage{
compmerge.list(l1, l2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l1}{a list}
  \item{l2}{a list}
}

\details{
\code{l1} and \code{l2} need to have same number of components.
}

\value{
Returns the merged list.
}

\author{
Sandipan Pramanik, Valen E. Johnson and Anirban Bhattacharya
}

\examples{
list1 = list("a"=1, "b"=3, "c"=5)
list2 = list("a"=2, "b"=4, "c"=6)

compmerge.list(list1, list2)
}

