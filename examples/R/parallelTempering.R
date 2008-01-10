
### $Id: parallelTempering.R,v 1.2 2008/01/11 02:31:14 goswami Exp $

if (!('package:EMC' %in% search( )))
    library(EMC, lib.loc='~/Rlib')

source('../../tests/funcGenerators.R')

pdf(file = 'parallelTempering.pdf')

samplerObj <-
    with(VShapedFuncGenerator( ),
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

dev.off( )
