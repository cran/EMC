
### $Id: placeTempers3.R,v 1.2 2008/01/11 02:31:14 goswami Exp $

if (!('package:EMC' %in% search( )))
    library(EMC, lib.loc='~/Rlib')

source('../../tests/funcGenerators.R')

pdf(file = 'placeTempers3.pdf')

placeTempersObj <-
    with(uniModeFuncGenerator( ),
         placeTempers(nIters            = 10000,
                      acceptRatioLimits = c(0.5, 0.6),
                      ladderLenMax      = 50,
                      startingVals      = c(0, 0),
                      logTarDensFunc    = logTarDensFunc,
                      MHPropNewFunc     = MHPropNewFunc,
                      temperLimits      = c(1, 3),
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
              xlim = c(-7, 7),
              ylim = c(-7, 7),
              pch  = '.',
              ask  = FALSE,
              main = as.expression(main),
              xlab = as.expression(substitute(x[xii], list(xii = 1))),
              ylab = as.expression(substitute(x[xii], list(xii = 2))))   
     }
 })


dev.off( )
