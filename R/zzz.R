.First.lib <- function(lib, pkg) {
    if(!require(exactRankTests))
        warning("Could not load package exactRankTests")  
    if(!require(mvtnorm))
        warning("Could not load package mvtnorm")
}