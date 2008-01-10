
### $Id: zzz.R,v 1.9 2008/01/11 02:31:13 goswami Exp $

.onLoad <-
    function (libname, pkgname)
{
    this.year <- substr(as.character(Sys.Date( )), 1, 4)    
    cat('##\n',
        '## Evolutionary Monte Carlo Package (EMC)\n',
        '##\n',
        '## Functionality: random walk Metropolis, general Metropolis-Hastings\n',
        '## parallel tempering, target oriented EMC (TOEMC), temperature ladder\n',
        '## contruction and placement\n',
        '##\n',
        '## Use: "help(package = EMC)" at the R prompt for more info\n',
        '##\n',
        '## Copyright (C) 2006-', this.year, ' Gopi Goswami\n',
        '##\n',
        '## Created and maintained by: Gopi Goswami <goswami@stat.harvard.edu>\n',
        '##\n', sep = '')        

    require(mvtnorm, quietly = TRUE)
    library.dynam(pkgname, pkgname, lib.loc = libname)
}






