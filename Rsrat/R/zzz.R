
.onLoad <- function(lib, pkg){
    packageStartupMessage( "Rsrat: A Software Reliability Assessment Tool" )
    library.dynam("Rsrat", pkg, lib)
}
