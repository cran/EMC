
%%% $Id: MetropolisHastings.Rd,v 1.18 2008/01/11 04:12:41 goswami Exp $

\name{MetropolisHastings}

\alias{MetropolisHastings}

\title{The Metropolis-Hastings algorithm}

\description{
  Given a target density function and an asymmetric proposal
  distribution, this function produces samples from the target using the
  Metropolis Hastings algorithm.
}

\usage{
MetropolisHastings(nIters,
                   startingVal,         
                   logTarDensFunc,      
                   propNewFunc,         
                   logPropDensFunc,     
                   MHBlocks      = NULL,
                   MHBlockNTimes = NULL,
                   saveFitness   = FALSE, 
                   verboseLevel  = 0,
                   \dots)    
}

\arguments{
  Below \code{sampDim} refers to the dimension of the sample space.

  \item{nIters}{\code{integer} \eqn{>} 0.}
  
  \item{startingVal}{\code{double} vector of length \code{sampDim}.}
  
  \item{logTarDensFunc}{\code{function} of two arguments
    \code{(draw, \dots)}
    that returns the target density evaluated in the log scale.}

  \item{propNewFunc}{\code{function} of three arguments
    \code{(block, currentDraw, \dots)}
    that returns new Metropolis-Hastings proposals. \emph{See details
      below on the argument block.}}

  \item{logPropDensFunc}{\code{function} of four arguments
    \code{(block, currentDraw, proposalDraw, \dots)} that returns the
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

  \item{verboseLevel}{\code{integer}, a value \eqn{\ge}{>=} 2 produces a
    lot of output.}

  \item{\dots}{optional arguments to be passed to \code{logTarDensFunc},
    \code{propNewFunc} and \code{logPropDensFunc}.}
}

\details{
  \describe{
    \item{\code{propNewFunc} and \code{logPropDensFunc}}{\cr
      The \code{propNewFunc} and the \code{logPropDensFunc} are called
      multiple times by varying the \code{block} argument over 
      \code{1:length(MHBlocks)}, so these functions should know how to
      generate a proposal from the \code{currentDraw} or to evaluate
      the proposal density depending on which block was passed as the
      argument. \emph{See the example section for sample code.}}

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
}

\value{
  This function returns a list with the following components:
  \item{draws}{\code{matrix} of dimension \code{nIters} \eqn{\times}{x}
      \code{sampDim}, if \code{saveFitness = FALSE}. If
      \code{saveFitness = TRUE}, then the returned matrix is of
      dimension \code{nIters} \eqn{\times}{x} \code{(sampDim + 1)},
      where the fitness values appear in its last column.}

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
  }
}

\author{Gopi Goswami \email{goswami@stat.harvard.edu}}

\references{
  \cite{Jun S. Liu (2001). \emph{Monte Carlo strategies for
      scientific computing.} Springer.}
}

\seealso{\code{\link{randomWalkMetropolis}},
  \code{\link{parallelTempering}}, \code{\link{evolMonteCarlo}}}

\examples{
## The Cigar-shaped distribution with independent t-proposal
CigarShapedFuncGenerator2 <-
    function (seed = 13579)
{
    set.seed(seed)    
    dd     <- 2
    ARDisp <-
        function (rho)
        {
            tmp <- rep(1, dd)
            diag((1 - rho) * tmp) + rho * tmp \%*\% t(tmp)
        }

    means          <- c(1, 1)
    disp           <- ARDisp(-0.95)
    logTarDensFunc <-
        function (draw, ...)
            dmvnorm(draw, means, disp, log = TRUE)

    tDF         <- 3
    propNewFunc <-
        function (block, currentDraw, ...)
        {
            proposalDraw        <- currentDraw
            proposalDraw[block] <- rt(1, tDF)
            proposalDraw
        }            
    
    logPropDensFunc <-
        function (block, currentDraw, proposalDraw, ...)
            dt(proposalDraw[block], tDF, log = TRUE)

    list(logTarDensFunc  = logTarDensFunc,
         propNewFunc     = propNewFunc,
         logPropDensFunc = logPropDensFunc)
}


samplerObj <-
    with(CigarShapedFuncGenerator2( ),
         MetropolisHastings(nIters          = 5000,
                            startingVal     = c(0, 0),
                            logTarDensFunc  = logTarDensFunc,
                            propNewFunc     = propNewFunc,
                            logPropDensFunc = logPropDensFunc,
                            verboseLevel    = 2))
print(samplerObj)
print(names(samplerObj))
with(samplerObj,
 {
     print(detailedAcceptRatios)
     print(dim(draws))
     plot(draws,
          xlim = c(-3, 5),
          ylim = c(-3, 4),
          pch  = '.',
          ask  = FALSE,
          main = as.expression(paste('# draws:', nIters)),
          xlab = as.expression(substitute(x[xii], list(xii = 1))),
          ylab = as.expression(substitute(x[xii], list(xii = 2))))    
 })
}

\keyword{methods}
