
%%% $Id: placeTempers.Rd,v 1.17 2008/02/04 17:01:38 goswami Exp $

\name{placeTempers}

\alias{placeTempers}

\title{Place the intermediate temperatures between the temperature limits}

\description{
  Multiple MCMC chains based algorithms (e.g., parallel tempering,
  evolutionary Monte Carlo) need a temperature ladder. This function
  places the intermediate temperatures between the minimum and the
  maximum temperature for the ladder.
} 

\usage{
placeTempers(nIters,                           
             acceptRatioLimits,                      
             ladderLenMax,                           
             startingVals,                           
             logTarDensFunc,                         
             MHPropNewFunc,                          
             logMHPropDensFunc = NULL,               
             temperLadder      = NULL,               
             temperLimits      = NULL,               
             ladderLen         = 15,                 
             scheme            = 'exponential',      
             schemeParam       = 1.5,                
             guideMe           = TRUE,               
             levelsSaveSampFor = NULL,               
             saveFitness       = FALSE,              
             verboseLevel      = 0,
             \dots)                  
}

\arguments{
  Below \code{sampDim} refers to the dimension of the sample space,
  \code{temperLadderLen} refers to the length of the temperature ladder,
  and \code{levelsSaveSampForLen} refers to the length of
  \code{levelsSaveSampFor}. Note, this function calls
  \code{\link{evolMonteCarlo}}, so some of the arguments below have the
  same name and meaning as the corresponding ones for
  \code{\link{evolMonteCarlo}}. \emph{See details below for explanation
  on the following arguments}.

  \item{nIters}{\code{integer} \eqn{>} 0.}

  \item{acceptRatioLimits}{\code{double} vector of \emph{two
      probabilities}.}

  \item{ladderLenMax}{\code{integer} \eqn{>} 0.}

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
  This function is based on the temperature placement method introduced
  in section 4.2 of Goswami and Liu (2007).

  \describe{ \item{\code{acceptRatioLimits}}{This is a range for the
    estimated acceptance ratios for the random exchange move for the
    consecutive temperature levels of the final ladder. It is
    recommended that specified range is between 0.3 and 0.6.}
    
    \item{\code{ladderLenMax}}{It is preferred that one specifies
      \code{acceptRatioLimits} for constructing the final temperature
      ladder. However, If one has some computational limitations then
      one could also specify \code{ladderLenMax} which will limit the
      length of the final temperature ladder produced. This also serves
      as an upper bound on the number of temperature levels while
      placing the intermediate temperatures using the
      \code{acceptRatioLimits}.}

    \item{\code{temperLadder}}{This is the temperature ladder needed for
      the second stage preliminary run. One can either specify a
      temperature ladder via \code{temperLadder} or specify
      \code{temperLimits}, \code{ladderLen}, \code{scheme} and
      \code{schemeParam}. For details on the later set of parameters,
      see below. Note, \code{temperLadder} overrides
      \code{temperLimits}, \code{ladderLen}, \code{scheme} and
      \code{schemeParam}.}
    
    \item{\code{temperLimits}}{\code{temperLimits = c(lowerLimit,
      upperLimit)} is a two-tuple of positive numbers, where the
      \code{lowerLimit} is usually 1 and \code{upperLimit} is a number
      in [100, 1000]. If stochastic optimization (via sampling) is the
      goal, then \code{lowerLimit} is taken to be in [0, 1]. Often the
      \code{upperLimit} is the maximum temperature as suggested by
      \code{\link{findMaxTemper}}.}

    \item{\code{ladderLen}, \code{scheme} and \code{schemeParam}}{These
      three parameters are required (along with \code{temperLimits}) if
      \code{temperLadder} \emph{is not} provided. We recommend taking
      \code{ladderLen} in [15, 30]. The allowed choices for
      \code{scheme} and \code{schemeParam} are:}

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
    
    \item{\code{guideMe}}{If \code{guideMe = TRUE}, then the function
      suggests different modifications to alter the setting towards a
      re-run, in case there are problems with the underlying MCMC run.}

    \item{\code{levelsSaveSampFor}}{This is passed to
      \code{evolMonteCarlo} for the underlying MCMC run.}

    }
}

\value{
  This function returns a list with the following components:

  \item{finalLadder}{the final temperature ladder found by placing the
    intermediate temperatures to be used in \code{parallelTempering} or
    \code{evolMonteCarlo}.}

  \item{temperLadder}{the temperature ladder used for the underlying
    MCMC run.}

  \item{acceptRatiosEst}{the estimated acceptance ratios for the random
      exchange move for the consecutive temperature levels of
      \code{temperLadder}.}

  \item{CVSqWeights}{this is the square of the coefficient of variation
    of the weights of the importance sampling estimators used to
    estimate the acceptance ratios, namely, \code{estAcceptRatios}.}
  
  \item{temperLimits}{the sorted \code{temperLimits} argument.}

  \item{acceptRatioLimits}{the sorted \code{acceptRatioLimits} argument.}

  \item{nIters}{the post burn-in \code{nIters}.}

  \item{levelsSaveSampFor}{the \code{levelsSaveSampFor} argument.}

  \item{draws}{\code{array} of dimension \code{nIters} \eqn{\times}{x}
      \code{sampDim} \eqn{\times}{x} \code{levelsSaveSampForLen}, if
      \code{saveFitness = FALSE}. If \code{saveFitness = TRUE}, then the
      returned array is of dimension \code{nIters} \eqn{\times}{x}
      \code{(sampDim + 1)} \eqn{\times}{x} \code{levelsSaveSampForLen};
      i.e., each of the \code{levelsSaveSampForLen} matrices contain the
      fitness values in their last column.}

  \item{startingVals}{the \code{startingVals} argument.}

  \item{time}{the time taken by the run.}

}

\note{
  The effect of leaving the default value \code{NULL} for some of the
  arguments above are as follows:

  \tabular{rl}{
    \code{logMHPropDensFunc}
    \tab the proposal density \code{MHPropNewFunc} is deemed symmetric.
    \cr

    \code{temperLadder}
    \tab valid \code{temperLimits}, \code{ladderLen}, \code{scheme} and
    \code{schemeParam}
    \cr

    \tab
    are provided, which are used to construct the \code{temperLadder}.
    \cr
    
    \code{temperLimits}
    \tab a valid \code{temperLadder} is provided.
    \cr

    \code{levelsSaveSampFor}
    \tab \code{temperLadderLen}.

  }
}

\author{Gopi Goswami \email{goswami@stat.harvard.edu}}

\references{
  \cite{Gopi Goswami and Jun S. Liu (2007). On learning strategies for
    evolutionary Monte Carlo. \emph{Statistics and Computing
      17:1:23-38}.}
}

\seealso{\code{\link{findMaxTemper}}, \code{\link{parallelTempering}},
  \code{\link{evolMonteCarlo}}}

\examples{
placeTempersObj <-
    with(VShapedFuncGenerator(-13579),
         placeTempers(nIters            = 10000,
                      acceptRatioLimits = c(0.5, 0.6),
                      ladderLenMax      = 50,
                      startingVals      = c(0, 0),
                      logTarDensFunc    = logTarDensFunc,
                      MHPropNewFunc     = MHPropNewFunc,
                      temperLimits      = c(1, 5),
                      ladderLen         = 10,
                      levelsSaveSampFor = seq_len(10),
                      verboseLevel      = 1))
print(placeTempersObj)
print(names(placeTempersObj))
with(placeTempersObj,
 {
     par(mfcol = c(3, 3))
     for (ii in seq_along(levelsSaveSampFor)) {
         main <- paste('temper:', round(temperLadder[levelsSaveSampFor[ii]], 3))
         plot(draws[ , , ii],
              xlim = c(-4, 20),
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
