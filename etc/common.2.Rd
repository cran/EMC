
\name{common.2}

\alias{common.2}

\title{The common stuff for parallelTempering, evolMonteCarlo}

\arguments{
  Below \code{sampDim} refers to the dimension of the sample space,
  \code{temperLadderLen} refers to the length of the temperature ladder,
  and \code{levelsSaveSampForLen} refers to the length of 
  \code{levelsSaveSampFor}. \emph{See details below for explanation on
    the following arguments}.
  
  \item{nIters}{\code{integer} \eqn{>} 0.}
  
  \item{temperLadder}{\code{double} vector with all \emph{positive}
    entries, in \emph{decreasing} order.}
  
  \item{startingVals}{\code{double} matrix of dimension
    \code{temperLadderLen} \eqn{\times}{x} \code{sampDim} or vector of
    length \code{sampDim}, in which case the same starting values are
    used for every temperature level.}
  
  \item{logTarDensFunc}{\code{function} of two arguments
    \code{(draw, ...)}
    that returns the target density evaluted in the log scale.}
  
  \item{MHPropNewFunc}{\code{function} of four arguments
    \code{(temperature, block, currentDraw, ...)}
    that returns new Metropolis-Hastings proposals. \emph{See details
      below on the argument block.}}

  \item{logMHPropDensFunc}{\code{function} of five arguments
    \code{(temperature, block, currentDraw, proposalDraw, ...)}
    that returns the proposal density evaluated in the log
    scale. \emph{See details below on the argument block.}}
  
  \item{MHBlocks}{\code{list} of integer vectors giving dimensions to be
    blocked together for sampling. It defaults to
    \code{as.list(1:sampDim)}, i.e., each dimension is treated as a
    block on its own. \emph{See details below for an example.}}

  \item{MHBlockNTimes}{\code{integer} vector of number of times each block
    given by \code{MHBlocks} should be sampled in each iteration. It
    defaults to \code{rep(1, length(MHBlocks))}. \emph{See details below
      for an example.}}

  \item{moveProbsList}{named \code{list} of probabilities adding upto 1.}

  \item{moveNTimesList}{named \code{list} of integers \eqn{\ge}{>=} 0.}

  \item{levelsSaveSampFor}{\code{integer} vector with \emph{positive}
    entries.}

  \item{saveFitness}{\code{logical}.}

  \item{verboseLevel}{\code{integer}, a value \eqn{\ge}{>=} 2 produces a
    lot of output.}
}

\details{
  \item{\code{MHPropNewFunc} and \code{logMHPropDensFunc}}{\cr
    The \code{MHPropNewFunc} and the \code{logMHPropDensFunc} are called
    multiple times by varying the \code{block} argument over 
    \code{1:length(MHBlocks)}, so these functions should know how to
    generate a proposal from the \code{currentDraw} or to evaluate
    the proposal density depending on which block of was passed as the
    argument. \emph{See the example section for sample code.}}

  \item{\code{MHBlocks} and \code{MHBlockNTimes}}{\cr
    Blocking is an important and useful tool in MCMC that helps
    speed up sampling and hence mixing. Example: Let \code{sampDim =
      6}. Let we want to sample dimensions 1, 2, 4 as one block,
    dimensions 3 and 5 as another and treat dimension 6 as the third
    block. Suppose we want to sample the three blocks mentioned above 1,
    5 and 10 times in each iteration, respectively. Then we could set
    \code{MHBlocks = list(c(1, 2, 4), c(3, 5), 6)} and
    \code{MHBlockNTimes = c(1, 5, 10)}}.
}

\value{
  This function returns a list with the following components:

  \item{draws}{\code{array} of dimension \code{nIters} \eqn{\times}{x}
      \code{sampDim} \eqn{\times}{x} \code{levelsSaveSampForLen}, if
      \code{saveFitness = FALSE}. If \code{saveFitness = TRUE}, then the
      returned array is of dimension \code{nIters} \eqn{\times}{x}
      \code{(sampDim + 1)} \eqn{\times}{x} \code{levelsSaveSampForLen};
      i.e., each of the \code{levelsSaveSampForLen} matrices contain the
      fitness values in their last column.}

  \item{acceptRatios}{\code{matrix} of the acceptance rates for various
    moves used.}

  \item{detailedAcceptRatios}{\code{list} of matrices with detailed
    summary of the acceptance rates for various moves used.}

  \item{nIters}{the \code{nIters} argument.}

  \item{temperLadder}{the \code{temperLadder} argument.}

  \item{startingVals}{the \code{startingVals} argument.}

  \item{moveProbsList}{the \code{moveProbsList} argument.}

  \item{moveNTimesList}{the \code{moveNTimesList} argument.}  

  \item{time}{the time taken by the run.}
}

\note{
  The effect of leaving the default value \code{NULL} for some of the
  arguments above are as follows:

  \tabular{rl}{
    \code{logMHPropDensFunc}
    \tab the proposal density \code{MHPropNewFunc} is deemed symmetric.
    \cr

    \code{MHBlocks}
    \tab \code{as.list(1:sampDim)}.
    \cr

    \code{MHBlockNTimes}
    \tab \code{rep(1, length(MHBlocks))}.
    \cr

    \code{moveProbesList}
    \tab \code{list(MH = 0.4, RC = 0.3, SC = 0.3)}.
    \cr

    \code{moveNTimesList}
    \tab \code{list(MH = 1, RC = mm, SC = mm, RE = nn)}, where 
    \cr

    \tab
    \code{mm <- floor(nn / 2)} and \code{nn <- temperLadderLen}.
    \cr

    \code{levelsSaveSampFor}
    \tab \code{temperLadderLen}.
  }
}

\author{Gopi Goswami \email{goswami@stat.harvard.edu}}

\references{
  \cite{Faming Liang and Wing H.Wong (2001). Real-Parameter Evolutionary
    Monte Carlo with Applications to Bayesian Mixture
    Models. \emph{Journal of the American Statistical Association
      96:653-666}.}
}



