
%%% $Id: parallelTempering.Rd,v 1.21 2008/02/04 17:01:38 goswami Exp $

\name{parallelTempering}

\alias{parallelTempering}

\title{The parallel Tempering algorithm}

\description{
  Given a multi-modal and multi-dimensional target density function, a
  (possibly asymmetric) proposal distribution and a temperature ladder,
  this function produces samples from the target using the parallel
  tempering algorithm.
}

\usage{
parallelTempering(nIters,                         
                  temperLadder,                   
                  startingVals,                   
                  logTarDensFunc,                 
                  MHPropNewFunc,                  
                  logMHPropDensFunc = NULL,       
                  MHBlocks          = NULL,
                  MHBlockNTimes     = NULL,
                  moveProbsList     = NULL,           
                  moveNTimesList    = NULL,          
                  levelsSaveSampFor = NULL,       
                  saveFitness       = FALSE,            
                  verboseLevel      = 0,
                  \dots)
}

\arguments{
  Below \code{sampDim} refers to the dimension of the sample space,
  \code{temperLadderLen} refers to the length of the temperature ladder,
  and \code{levelsSaveSampForLen} refers to the length of the
  \code{levelsSaveSampFor}.
  
  \item{nIters}{\code{integer} \eqn{>} 0.}
  
  \item{temperLadder}{\code{double} vector with all \emph{positive}
    entries, in \emph{decreasing} order.}
  
  \item{startingVals}{\code{double} matrix of dimension
    \code{temperLadderLen} \eqn{\times}{x} \code{sampDim} or vector of
    length \code{sampDim}, in which case the same starting values are
    used for every temperature level.}
  
  \item{logTarDensFunc}{\code{function} of two arguments
    \code{(draw, \dots)}
    that returns the target density evaluated in the log scale.}
  
  \item{MHPropNewFunc}{\code{function} of four arguments
    \code{(temperature, block, currentDraw, \dots)}
    that returns new Metropolis-Hastings proposals. \emph{See details
      below on the argument block.}}

  \item{logMHPropDensFunc}{\code{function} of five arguments
    \code{(temperature, block, currentDraw, proposalDraw, \dots)}
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

  \item{\dots}{optional arguments to be passed to \code{logTarDensFunc},
    \code{MHPropNewFunc} and \code{logMHPropDensFunc}.}
  
}

\details{
  \describe{

    \item{\code{MHPropNewFunc} and \code{logMHPropDensFunc}}{The
      \code{MHPropNewFunc} and the \code{logMHPropDensFunc} are called
      multiple times by varying the \code{block} argument over
      \code{1:length(MHBlocks)}, so these functions should know how to
      generate a proposal from the \code{currentDraw} or to evaluate the
      proposal density depending on which block was passed as the
      argument. \emph{See the example section for sample code.}}

    \item{\code{MHBlocks} and \code{MHBlockNTimes}}{Blocking is an
      important and useful tool in MCMC that helps speed up sampling and
      hence mixing. Example: Let \code{sampDim = 6}. Let we want to
      sample dimensions 1, 2, 4 as one block, dimensions 3 and 5 as
      another and treat dimension 6 as the third block. Suppose we want
      to sample the three blocks mentioned above 1, 5 and 10 times in
      each iteration, respectively. Then we could set \code{MHBlocks =
      list(c(1, 2, 4), c(3, 5), 6)} and \code{MHBlockNTimes = c(1, 5,
      10)}}.

    \item{The parallel tempering algorithm}{The parallel tempering (PT;
      Liang and Wong, 2001) algorithm is composed of the following
      moves:}    

    \tabular{rl}{
      \acronym{MH} \tab Metropolis-Hastings or mutation \cr
      \acronym{RE} \tab (random) exchange 
    }

    The current function could be used to run the PT algorithm by
    specifying what moves to employ using the following variables.

    \describe{

      \item{\code{moveProbsList} and \code{moveNTimesList}}{The allowed
        names for components of \code{moveProbsList} and
        \code{moveNTimesList} come from the abbreviated names of the
        moves above. For example, the following specifications are
        valid:}

      \preformatted{
moveProbsList = list(MH = 0.4,
                     RE = 0.6)
       }

       \preformatted{
moveNTimesList = list(MH = 1,
                      RE = temperLadderLen)           
        }
    }    

    \item{\code{levelsSaveSampFor}}{By default, samples are saved and
      returned for temperature level \code{temperLadderLen}. The
      \code{levelsSaveSampFor} could be used to save samples from other
      temperature levels as well (e.g., \code{levelsSaveSampFor =
      1:temperLadderLen} saves samples from all levels).}

    \item{\code{saveFitness}}{The term \emph{fitness} refers to the
      function \eqn{H(x)}, where the target density of interest is given
      by:

      \deqn{g(x) \propto \exp[ -H(x) / \tau_{min} ]}

      \eqn{H(x)} is also known as the \emph{energy} function. By default,
      the fitness values are not saved, but one can do so by setting
      \code{saveFitness = TRUE}.}
    
  }
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

  \item{levelsSaveSampFor}{the \code{levelsSaveSampFor} argument.}

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

    \code{moveProbsList}
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

\seealso{\code{\link{evolMonteCarlo}}}

\examples{
samplerObj <-
    with(VShapedFuncGenerator(-13579),
         parallelTempering(nIters            = 2000,
                           temperLadder      = c(15, 6, 2, 1),
                           startingVals      = c(0, 0),
                           logTarDensFunc    = logTarDensFunc,
                           MHPropNewFunc     = MHPropNewFunc,
                           levelsSaveSampFor = seq_len(4), 
                           verboseLevel      = 1))
print(samplerObj)
print(names(samplerObj))
with(samplerObj,
 {
     print(detailedAcceptRatios)
     print(dim(draws))
     par(mfcol = c(2, 2))
     for (ii in seq_along(levelsSaveSampFor)) {
         main <- paste('temper:', round(temperLadder[levelsSaveSampFor[ii]], 3))
         plot(draws[ , , ii],
              xlim = c(-5, 20),
              ylim = c(-8, 8),
              pch  = '.',
              ask  = FALSE,
              main = as.expression(main),
              xlab = as.expression(substitute(x[xii], list(xii = 1))),
              ylab = as.expression(substitute(x[xii], list(xii = 2))))
     }
 })
}

\keyword{methods}
