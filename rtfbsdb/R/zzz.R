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
	exit_code <- system("starch --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* The starch command in bedops package doesn't work normally.")

	exit_code <- system("starchcat --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* The starchcat command in bedops package doesn't work normally.")

	exit_code <- system("sort-bed --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* The sort-bed command in bedops package doesn't work normally.")

	exit_code <- system("awk --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* The awk command in bedops package doesn't work normally.")

	exit_code <- system("zcat --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
		packageStartupMessage("* The zcat command in bedops package doesn't work normally.")

	exit_code <- system("twoBitInfo",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code == 127 )
		packageStartupMessage("* The twoBitInfo command in UCSC package doesn't work normally.")

	# dont use samtools any more since 4/2/2016
	# exit_linux <- system("samtools --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	# exit_osx <- system("samtools",  ignore.stderr=T, ignore.stdout=T, wait = T );
	# if( exit_linux != 0  && exit_osx != 0 )
	# 	packageStartupMessage("* The samtools command doesn't work normally (http://samtools.sourceforge.net/).")

	exit_linux <- system("bedtools --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	exit_osx <- system("bedtools",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_linux != 0  && exit_osx  != 0 )
		packageStartupMessage("* The bedtools command doesn't work normally (http://bedtools.readthedocs.org/).")
}