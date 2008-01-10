
### $Id: findMaxTemper1.R,v 1.2 2008/01/11 02:31:13 goswami Exp $

if (!('package:EMC' %in% search( )))
    library(EMC, lib.loc='~/Rlib')

source('../../tests/funcGenerators.R')

pdf(file = 'findMaxTemper1.pdf')

coord1        <- function (xx) { xx[1] }
ss            <- function (xx) { sum(xx) }
pp            <- function (xx) { prod(xx) }
statsFuncList <- list(coord1, ss, pp)
maxTemperObj  <-
    with(VShapedFuncGenerator( ),
         findMaxTemper(nIters            = 10000,
                       statsFuncList     = statsFuncList,
                       temperLadder      = c(20, 15, 10, 5, 1),
                       startingVals      = c(0, 0),
                       logTarDensFunc    = logTarDensFunc,
                       MHPropNewFunc     = MHPropNewFunc,
                       levelsSaveSampFor = seq_len(5),
                       doFullAnal        = TRUE,
                       verboseLevel      = 1))
print(maxTemperObj)
print(names(maxTemperObj))
with(maxTemperObj,
 {
     par(mfcol = c(3, 3))
     for (ii in seq_along(levelsSaveSampFor)) {
         main <- paste('temper:', round(temperLadder[levelsSaveSampFor[ii]], 3))
         plot(draws[ , , ii],
              xlim = c(-10, 25),
              ylim = c(-10, 10),
              pch  = '.',
              ask  = FALSE,
              main = as.expression(main),
              xlab = as.expression(substitute(x[xii], list(xii = 1))),
              ylab = as.expression(substitute(x[xii], list(xii = 2))))   
     }
 })

dev.off( )
