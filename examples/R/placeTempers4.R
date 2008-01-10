
### $Id: placeTempers4.R,v 1.2 2008/01/11 02:31:15 goswami Exp $

if (!('package:EMC' %in% search( )))
    library(EMC, lib.loc='~/Rlib')

source('../../tests/funcGenerators.R')

pdf(file = 'placeTempers4.pdf')

placeTempersObj <-
    with(twentyModeFuncGenerator( ),
         placeTempers(nIters            = 10000,
                      acceptRatioLimits = c(0.5, 0.6),
                      ladderLenMax      = 50,
                      startingVals      = c(0, 0),
                      logTarDensFunc    = logTarDensFunc,
                      MHPropNewFunc     = MHPropNewFunc,
                      temperLimits      = c(1, 15),
                      ladderLen         = 15,
                      levelsSaveSampFor = seq_len(15),
                      verboseLevel      = 2))
print(placeTempersObj)
print(names(placeTempersObj))
with(placeTempersObj,
 {
     par(mfcol = c(3, 3))
     for (ii in seq_along(levelsSaveSampFor)) {
         main <- paste('temper:', round(temperLadder[levelsSaveSampFor[ii]], 3))
         plot(draws[ , , ii],
              xlim = c(-3, 15),
              ylim = c(-3, 15),
              pch  = '.',
              ask  = FALSE,
              main = as.expression(main),
              xlab = as.expression(substitute(x[xii], list(xii = 1))),
              ylab = as.expression(substitute(x[xii], list(xii = 2))))   
     }
 })


dev.off( )
