
### $Id: groupTest1.R,v 1.2 2008/01/11 02:31:14 goswami Exp $

if (!('package:EMC' %in% search( )))
    library(EMC, lib.loc='~/Rlib')

source('../../tests/funcGenerators.R')

testEvolMonteCarlo <-
    function (controlsObj)
{
    samplerObj <-
        with(controlsObj,
         {
             cat(rep('=', 80), '\n',
                 'BEGIN: Testing ', exampleName, '\n\n',
                 sep = '')
             evolMonteCarlo(nIters            = nIters,
                            temperLadder      = temperLadder,
                            startingVals      = startingVals,
                            logTarDensFunc    = logTarDensFunc,
                            MHPropNewFunc     = MHPropNewFunc,
                            logMHPropDensFunc = logMHPropDensFunc,
                            moveProbsList     = moveProbsList,
                            moveNTimesList    = moveNTimesList, 
                            SCRWMNTimes       = SCRWMNTimes,
                            SCRWMPropSD       = SCRWMPropSD,
                            levelsSaveSampFor = levelsSaveSampFor,
                            verboseLevel      = verboseLevel)
         })

    print(samplerObj)
    print(names(samplerObj))
    print(samplerObj$acceptRatios)
    print(samplerObj$detailedAcceptRatios)
    print(dim(samplerObj$draws))
    print(samplerObj$time)
    temperLadder      <- samplerObj$temperLadder
    nIters            <- samplerObj$nIters
    levelsSaveSampFor <- samplerObj$levelsSaveSampFor
    
    pp       <- samplerObj$moveProbsList[samplerObj$moveProbsList > 0.0]
    tt       <- samplerObj$moveNTimesList[samplerObj$moveNTimesList > 0]
    plotFile <- paste(controlsObj$exampleName, '-',
                      paste(names(pp), pp, sep = '', collapse = '-'), '-',
                      paste(names(tt), tt, sep = '', collapse = '-'),
                      '.pdf', sep = '')
    pdf(plotFile)
    par(mfcol = c(2, 2))
    for (ii in seq_along(levelsSaveSampFor)) {
        temper <- round(temperLadder[levelsSaveSampFor[ii]], 3)
        plot(samplerObj$draws[ , , ii],
             xlim = controlsObj$xlim,
             ylim = controlsObj$ylim,
             pch  = '.',
             ask  = FALSE,
             main = paste('temper:', temper, ', # draws:', nIters),
             xlab = as.expression(substitute(x[xii], list(xii = 1))),
             ylab = as.expression(substitute(x[xii], list(xii = 2))))
    }
    dev.off( )
    cat('\nE N D: Testing', controlsObj$exampleName, '\n\n')    
}

testWhapedAll <-
    function (seed = 13579)
{
    WShapedObj <- WShapedFuncGenerator(seed)
    paramsList <- list(nIters            = 5000,
                       temperLadder      = c(15, 6, 2, 1),
                       logMHPropDensFunc = NULL,
                       startingVals      = c(0, 0),
                       SCRWMNTimes       = 1,
                       SCRWMPropSD       = 3.0,
                       levelsSaveSampFor = NULL,
                       verboseLevel      = 1,
                       xlim              = c(-20, 20),
                       ylim              = c(-8, 8))

    controlsObj <-
        c(WShapedObj,
          paramsList,
          list(exampleName            = 'RCRandomRandom',
               moveSelectionCodesList = list(RC = c('random', 'random')),
               moveProbsList          = list(MH = 0.3, RC = 0.7),
               moveNTimesList         = list(MH = 1, RC = 1))
          )
    testEvolMonteCarlo(controlsObj)

    controlsObj <-
        c(WShapedObj,
          paramsList,
          list(exampleName            = 'RCBestBest',
               moveSelectionCodesList = list(RC = c('best', 'best')),
               moveProbsList          = list(MH = 0.3, RC = 0.7),
               moveNTimesList         = list(MH = 1, RC = 1))
          )
    testEvolMonteCarlo(controlsObj)
    
}


testWhapedAll( )


