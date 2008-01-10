
### $Id: evolMonteCarlo1.R,v 1.2 2008/01/11 02:31:13 goswami Exp $

if (!('package:EMC' %in% search( )))
    library(EMC, lib.loc='~/Rlib')

source('../../tests/funcGenerators.R')

pdf(file = 'evolMonteCarlo1.pdf')

samplerObj <-
    with(VShapedFuncGenerator( ),
     {
         allMovesNTimesList <- list(MH = 1, RC = 2, SC = 2, RE = 4,
                                    BIRE = 2, BCE = 2, BSE = 2, CE = 2)
         evolMonteCarlo(nIters            = 2000,
                        temperLadder      = c(15, 6, 2, 1),
                        startingVals      = c(0, 0),
                        logTarDensFunc    = logTarDensFunc,
                        MHPropNewFunc     = MHPropNewFunc,
                        moveNTimesList    = allMovesNTimesList,
                        SCRWMNTimes       = 1,
                        SCRWMPropSD       = 3.0,
                        levelsSaveSampFor = seq_len(4),
                        verboseLevel      = 1)
     })
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
