
### $Id: MetropolisHastings.R,v 1.2 2008/01/11 02:31:13 goswami Exp $

if (!('package:EMC' %in% search( )))
    library(EMC, lib.loc='~/Rlib')

source('../../tests/funcGenerators.R')

pdf(file = 'MetropolisHastings.pdf')

samplerObj <-
    with(CigarShapedFuncGenerator2( ),
         MetropolisHastings(nIters          = 5000,
                            startingVal     = c(0, 0),
                            logTarDensFunc  = logTarDensFunc,
                            propNewFunc     = propNewFunc,
                            logPropDensFunc = logPropDensFunc,
                            verboseLevel    = 2))
print(samplerObj)
print(names(samplerObj))
with(samplerObj,
 {
     print(detailedAcceptRatios)
     print(dim(draws))
     plot(draws,
          xlim = c(-3, 5),
          ylim = c(-3, 4),
          pch  = '.',
          ask  = FALSE,
          main = as.expression(paste('# draws:', nIters)),
          xlab = as.expression(substitute(x[xii], list(xii = 1))),
          ylab = as.expression(substitute(x[xii], list(xii = 2))))    
 })

dev.off( )
