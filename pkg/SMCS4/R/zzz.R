
.onLoad <-
    function (libname, pkgname)
{
    library.dynam(pkgname, pkgname, lib.loc=libname)
}






