.onLoad <- function(lib, pkg) {
    if(!require(survival))
        warning("Could not load package survival")
}
