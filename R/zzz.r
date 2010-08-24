# first and last
.onLoad <- function(libname,pkgname)
{
require("utils",quietly=TRUE,warn.conflicts=FALSE)
require("stats",quietly=TRUE,warn.conflicts=FALSE)
require("foreign",quietly=TRUE,warn.conflicts=FALSE)
require("gam",quietly=TRUE,warn.conflicts=FALSE)
require("splines",quietly=TRUE,warn.conflicts=FALSE)
require("gplots",quietly=TRUE,warn.conflicts=FALSE)
ver <- packageDescription(pkgname,fields="Version")
cat("This is \'ares\' library",ver,"\n",sep=" ")
}

.Last.lib <- function(libpath)
{
# nothing so far
}
