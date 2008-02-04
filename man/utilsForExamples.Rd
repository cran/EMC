
%%% $Id: utilsForExamples.Rd,v 1.1 2008/02/04 18:06:12 goswami Exp $

\name{utilsForExamples}

\alias{CigarShapedFuncGenerator1}
\alias{CigarShapedFuncGenerator2}
\alias{VShapedFuncGenerator}
\alias{WShapedFuncGenerator}
\alias{uniModeFuncGenerator}
\alias{twentyModeFuncGenerator}

\title{The utility function(s) for examples}

\description{
  The utility function(s) that are used in the example sections of the
  exported functions in this package.
}

\usage{
CigarShapedFuncGenerator1(seed)
CigarShapedFuncGenerator2(seed)
VShapedFuncGenerator(seed)
WShapedFuncGenerator(seed)
uniModeFuncGenerator(seed)
twentyModeFuncGenerator(seed)
}

\arguments{
  \item{seed}{the seed for random number generation.}
}

\value{
  A list containing the objects to be used as arguments to the exported
  functions in the respective example sections of this package. 
}

\author{Gopi Goswami \email{goswami@stat.harvard.edu}}

\seealso{\code{\link{randomWalkMetropolis}},
  \code{\link{MetropolisHastings}}, \code{\link{parallelTempering}},
  \code{\link{evolMonteCarlo}}, \code{\link{findMaxTemper}},
  \code{\link{placeTempers}}}

\keyword{datagen}

