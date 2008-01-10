
\name{common.3}

\alias{common.3}

\title{The common stuff for findMaxTemper, placeTempers}

\arguments{
  Below \code{sampDim} refers to the dimension of the sample space,
  \code{temperLadderLen} refers to the length of the temperature ladder,
  and \code{levelsSaveSampForLen} refers to the length of 
  \code{levelsSaveSampFor}. Note, this function calls
  \code{evolMonteCarlo}, so some of the arguments below have the same
  name and meaning as the corresponding ones for
  \code{evolMonteCarlo}. \emph{See details below for explanation on
    the following arguments}.

  \item{nIters}{\code{integer} \eqn{>} 0.}

  \item{startingVals}{\code{double} matrix of dimension
    \code{temperLadderLen} \eqn{\times}{x} \code{sampDim} or vector of
    length \code{sampDim}, in which case the same starting values are
    used for every temperature level.}

  \item{logTarDensFunc}{\code{function} of two arguments
    \code{(draw, \dots)} that returns the target density evaluated in
    the log scale.}
  
  \item{MHPropNewFunc}{\code{function} of four arguments
    \code{(temperature, block, currentDraw, \dots)}
    that returns new Metropolis-Hastings proposals. \emph{See details
      below on the argument block.}}

  \item{logMHPropDensFunc}{\code{function} of five arguments
    \code{(temperature, block, currentDraw, proposalDraw, \dots)}
    that returns the proposal density evaluated in the log
    scale. \emph{See details below on the argument block.}}
  
  \item{temperLadder}{\code{double} vector with all \emph{positive}
    entries, in \emph{decreasing} order.}

  \item{temperLimits}{\code{double} vector with \emph{two positive}
    entries.}

  \item{ladderLen}{\code{integer} \eqn{>} 0.}

  \item{scheme}{\code{character}.}

  \item{schemeParam}{\code{double} \eqn{>} 0.}

  \item{guideMe}{\code{logical}.}

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
    \item{\code{temperLimits}}{\cr
      \code{temperLimits = c(lowerLimit, upperLimit)} is a two-tuple
      of positive numbers, where the \code{lowerLimit} is usually 1
      and \code{upperLimit} is a number in [100, 1000]. If stochastic
      optimization (via sampling) is the goal, then \code{lowerLimit} is
      taken to be in [0, 1]. Often the \code{upperLimit} is the maximum
      temperature as suggested by \code{\link{findMaxTemper}}.}

    \item{\code{ladderLen}, \code{scheme} and \code{schemeParam}}{\cr
      These three parameters are required (along with
      \code{temperLimits}) if \code{temperLadder} \emph{is not}
      provided. We recommend taking \code{ladderLen} in [15, 30]. The
      allowed choices for \code{scheme} and \code{schemeParam} are:}

    \tabular{rl}{
      \code{scheme} \tab \code{schemeParam} \cr
      ======== \tab ============= \cr
      linear \tab NA \cr
      log \tab NA \cr
      geometric \tab NA \cr
      mult-power \tab NA \cr
      add-power \tab \eqn{\ge}{>=} 0 \cr
      reciprocal \tab NA \cr
      exponential \tab \eqn{\ge}{>=} 0 \cr
      tangent \tab \eqn{\ge}{>=} 0      
    }

    We recommended using \code{scheme = 'exponential'} and
    \code{schemeParam} in [1.5, 2].
    
    \item{\code{guideMe}}{\cr
      If \code{guideMe = TRUE}, then the function suggests different
      modifications to alter the setting towards a re-run, in case there
      are problems with the underlying MCMC run.}
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

  \item{temperLadder}{the temperature ladder used for the underlying
    MCMC run.}

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



