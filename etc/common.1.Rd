
\name{common.1}

\alias{common.1}

\title{The Common stuff for randomWalkMetropolis, MetropolisHastings}

\arguments{
  Below \code{sampDim} refers to the dimension of the sample space.

  \item{nIters}{\code{integer} \eqn{>} 0.}
  
  \item{startingVal}{\code{double} vector of length \code{sampDim}.}
  
  \item{logTarDensFunc}{\code{function} of two arguments
    \code{(draw, ...)}
    that returns the target density evaluated in the log scale.}

  \item{propNewFunc}{\code{function} of three arguments
    \code{(block, currentDraw, ...)}
    that returns new Metropolis-Hastings proposals. \emph{See details
      below on the argument block.}}

  \item{logPropDensFunc}{\code{function} of four arguments
    \code{(block, currentDraw, proposalDraw, ...)} that returns the
    proposal density evaluated in the log scale. \emph{See details below
      on the argument block.}}
  
  \item{MHBlocks}{\code{list} of integer vectors giving dimensions to be
    blocked together for sampling. It defaults to
    \code{as.list(1:sampDim)}, i.e., each dimension is treated as a
    block on its own. \emph{See details below for an example.}}

  \item{MHBlockNTimes}{\code{integer} vector of number of times each block
    given by \code{MHBlocks} should be sampled in each iteration. It
    defaults to \code{rep(1, length(MHBlocks))}. \emph{See details below
      for an example.}}

  \item{saveFitness}{\code{logical} indicating whether fitness values
    should be saved. \emph{See details below.}}

  \item{verboseLevel}{\code{integer} in [0, 3], a value \eqn{\ge}{>=} 2
    produces a lot of output.}
}


\details{
  \item{\code{MHBlocks} and \code{MHBlockNTimes}}{\cr
    Blocking is an important and useful tool in MCMC that helps
    speed up sampling and hence mixing. Example: Let \code{sampDim =
      6}. Let we want to sample dimensions 1, 2, 4 as one block,
    dimensions 3 and 5 as another and treat dimension 6 as the third
    block. Suppose we want to sample the three blocks mentioned above 1,
    5 and 10 times in each iteration, respectively. Then we could set
    \code{MHBlocks = list(c(1, 2, 4), c(3, 5), 6)} and
    \code{MHBlockNTimes = c(1, 5, 10)}}    

  \item{\code{saveFitness}}{\cr
    The term \emph{fitness} refers to the negative of the
    \code{logTarDensFunc} values. By default, the fitness values are
    not saved, but one can do so by setting \code{saveFitness =
      TRUE}.}
}

\value{
  This function returns a list with the following components:

  \item{draws}{\code{matrix} of dimension \code{nIters} \eqn{\times}{x}
      \code{sampDim}, if \code{saveFitness = FALSE}. If \code{saveFitness =
      TRUE}, then the returned matrix is of dimension \code{nIters}
      \eqn{\times}{x} \code{(sampDim + 1)}, where the fitness values
      appear in its last column.}
  
  \item{acceptRatios}{\code{matrix} of the acceptance rates.}
  
  \item{detailedAcceptRatios}{\code{matrix} with detailed summary of the
    acceptance rates.}
  
  \item{nIters}{the \code{nIters} argument.}
  
  \item{startingVal}{the \code{startingVal} argument.}
  
  \item{time}{the time taken by the run.}
}

\note{
  The effect of leaving the default value \code{NULL} for some of the
  arguments above are as follows:

  \tabular{rl}{
    \code{MHBlocks}
    \tab \code{as.list(1:sampDim)}.
    \cr

    \code{MHBlockNTimes}
    \tab \code{rep(1, length(MHBlocks))}.
    \cr
  }
}

\author{Gopi Goswami \email{goswami@stat.harvard.edu}}

\references{
  \cite{Jun S. Liu (2001). \emph{Monte Carlo strategies for
      scientific computing.} Springer.}
}


  
