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
	err_cmds<- c();
	if( exit_code != 0 )
	{
		#packageStartupMessage("* Can't find 'awk' comand which is required in some functions, please use Sys.which('awk') to locate it.")
		err_cmds <- c(err_cmds, "awk");
	}

	exit_code <- system("zcat --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
	{
		#packageStartupMessage("* Can't find 'zcat' comand which is required in some functions, please use Sys.which('zcat') to locate it.")
		err_cmds <- c(err_cmds, "zcat");
	}

	exit_code <- system("starch --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
	{
		#packageStartupMessage("* Can't find 'starch' comand which is required in some functions, please use Sys.which('starch') to locate it or install it.(http://bedops.readthedocs.org/en/latest/index.html)")
		err_cmds <- c(err_cmds, "starch");
	}

	exit_code <- system("starchcat --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
	{
		#packageStartupMessage("* Can't find 'starchcat' comand which is required in some functions, please use Sys.which('starchcat') to locate it or install it.(http://bedops.readthedocs.org/en/latest/index.html)")
		err_cmds <- c(err_cmds, "starchcat");
	}

	exit_code <- system("sort-bed --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code != 0 )
	{
		#packageStartupMessage("* Can't find 'sort-bed ' comand which is required in some functions, please use Sys.which('sort-bed ') to locate it or install it.(http://bedops.readthedocs.org/en/latest/index.html)")
		err_cmds <- c(err_cmds, "sort-bed");
	}

	exit_code <- system("twoBitInfo",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code == 127 )
	{
		#packageStartupMessage("* Can't find 'twoBitInfo' comand which is required in some functions, please use Sys.which('twoBitInfo') to locate it or install it.(http://hgdownload.cse.ucsc.edu/admin/exe/)")
		err_cmds <- c(err_cmds, "twoBitInfo");
	}

	exit_code <- system("twoBitToFa",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_code == 127 )
	{
		#packageStartupMessage("* Can't find 'twoBitToFa' comand which is required in some functions, please use Sys.which('twoBitToFa') to locate it or install it.(http://hgdownload.cse.ucsc.edu/admin/exe/)")
		err_cmds <- c(err_cmds, "twoBitToFa");
	}

	exit_linux <- system("bedtools --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	exit_osx <- system("bedtools",  ignore.stderr=T, ignore.stdout=T, wait = T );
	if( exit_linux != 0  && exit_osx  != 0 )
	{
		#packageStartupMessage("* Can't find 'bedtools' comand which is required in some functions, please use Sys.which('bedtools') to locate it or install it.(http://bedtools.readthedocs.org/)")
		err_cmds <- c(err_cmds, "bedtools");
	}

	# dont use samtools any more since 4/2/2016
	# exit_linux <- system("samtools --version",  ignore.stderr=T, ignore.stdout=T, wait = T );
	# exit_osx <- system("samtools",  ignore.stderr=T, ignore.stdout=T, wait = T );
	# if( exit_linux != 0  && exit_osx != 0 )
	# 	packageStartupMessage("* The samtools command doesn't work normally (http://samtools.sourceforge.net/).")

	if(length(err_cmds)>0)
	{
		packageStartupMessage("WARNING! The following dependencies were not found in your current environment.");
		packageStartupMessage("------");
		packageStartupMessage(paste(err_cmds, collapse=","));
		packageStartupMessage("------");
		packageStartupMessage("These dependencies are required by RTFBSDB to improve the processing speed.");
		packageStartupMessage("To troubleshoot these please do the following");
		packageStartupMessage("1. Make sure these commands are installed.");
		packageStartupMessage("2. check the $PATH variable to call Sys.getenv('PATH')");
		packageStartupMessage("   - or - check the dependency location by the command 'which', e.g. Sys.which('ls')");
		packageStartupMessage("3. If the command path can not be found in the default setting, plase call Sys.setenv() to set your command path.");
		packageStartupMessage("   e.g. Sys.setenv( PATH = paste(Sys.getenv('PATH'), '/your/command/path', sep=':') )" );
	}
}