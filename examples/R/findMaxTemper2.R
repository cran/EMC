
### $Id: findMaxTemper2.R,v 1.2 2008/01/11 02:31:13 goswami Exp $

if (!('package:EMC' %in% search( )))
    library(EMC, lib.loc='~/Rlib')

source('../../tests/funcGenerators.R')

pdf(file = 'findMaxTemper2.pdf')

coord1        <- function (xx) { xx[1] }
ss            <- function (xx) { sum(xx) }
pp            <- function (xx) { prod(xx) }
statsFuncList <- list(coord1, ss, pp)
maxTemperObj  <-
    with(WShapedFuncGenerator( ),
         findMaxTemper(nIters            = 10000,
                       statsFuncList     = statsFuncList,
                       startingVals      = c(0, 0),
                       logTarDensFunc    = logTarDensFunc,
                       MHPropNewFunc     = MHPropNewFunc,
                       temperLimits      = c(1, 100),
                       ladderLen         = 13,
                       levelsSaveSampFor = seq_len(13),
                       verboseLevel      = 2))
print(maxTemperObj)
print(names(maxTemperObj))
with(maxTemperObj,
 {
     par(mfcol = c(3, 3))
     for (ii in seq_along(levelsSaveSampFor)) {
         main <- paste('temper:', round(temperLadder[levelsSaveSampFor[ii]], 3))
         plot(draws[ , , ii],
              xlim = c(-20, 30),
              ylim = c(-20, 20),
              pch  = '.',
              ask  = FALSE,
              main = as.expression(main),
              xlab = as.expression(substitute(x[xii], list(xii = 1))),
              ylab = as.expression(substitute(x[xii], list(xii = 2))))   
     }
 })

dev.off( )
