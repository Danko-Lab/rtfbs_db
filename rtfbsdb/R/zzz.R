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

.onAttach<- function(libname, pkgName)
{
	exit_code <- system("awk --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* Can't find 'awk' comand which is required in some functions, please use Sys.which('awk') to locate it.")

	exit_code <- system("zcat --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* Can't find 'zcat' comand which is required in some functions, please use Sys.which('zcat') to locate it.")

	exit_code <- system("starch --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* Can't find 'starch' comand which is required in some functions, please use Sys.which('starch') to locate it or install it.(http://bedops.readthedocs.org/en/latest/index.html)")

	exit_code <- system("starchcat --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* Can't find 'starchcat' comand which is required in some functions, please use Sys.which('starchcat') to locate it or install it.(http://bedops.readthedocs.org/en/latest/index.html)")

	exit_code <- system("sort-bed --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* Can't find 'sort-bed ' comand which is required in some functions, please use Sys.which('sort-bed ') to locate it or install it.(http://bedops.readthedocs.org/en/latest/index.html)")

	exit_code <- system("twoBitInfo",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code == 127 )
		packageStartupMessage("* Can't find 'twoBitInfo' comand which is required in some functions, please use Sys.which('twoBitInfo') to locate it or install it.(http://hgdownload.cse.ucsc.edu/admin/exe/)")

	exit_linux <- system("bedtools --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	exit_osx <- system("bedtools",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_linux != 0  && exit_osx  != 0 )
		packageStartupMessage("* Can't find 'bedtools' comand which is required in some functions, please use Sys.which('bedtools') to locate it or install it.(http://bedtools.readthedocs.org/)")

	# dont use samtools any more since 4/2/2016
	# exit_linux <- system("samtools --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	# exit_osx <- system("samtools",  ignore.stderr=T, ignore.stdout=T, wait = T );
	# if( exit_linux != 0  && exit_osx != 0 )
	# 	packageStartupMessage("* The samtools command doesn't work normally (http://samtools.sourceforge.net/).")

}