.First.lib <- function(lib, pkg) {
    library.dynam("maxstat", pkg, lib)
    if(!require(exactRankTests))
        warning("Could not load package exactRankTests")  
    if(!require(mvtnorm))
        warning("Could not load package mvtnorm")
    if(!require(survival))
        warning("Could not load package survival")
}