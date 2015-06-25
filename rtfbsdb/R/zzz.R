#.First.lib <-
#    function(libname, pkgname, where)
#    library.dynam("atif", pkgname, libname)
#
#.Last.lib <-
#    function(libpath)
#    dyn.unload(file.path(libpath,
#                         "libs",
#                         paste("atif",
#                               .Platform$"dynlib.ext",
#                               sep = "")))
