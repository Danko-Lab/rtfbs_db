# Format:  pwm.matrix
#          transfac 
#          jaspar 
#          mscan
#          meme	

tfbs_importMotifs <- function(tfbs, format, filenames, motif_ids=NULL, skip.lines=0, pseudocount= -7, force_even= FALSE, ...)
{
	pwm.matrice <- list();
	pwm.filenames <- c();
	
	if( length(format)==1 && format=="pwm.matrix" )
	{
		## the paramater 'filenames' is a string vector. searching wildchars for the extension file "*.pwm".
		## by default, the filename will be motif ID.
		
		loc.wildchar <- regexpr("\\*",filenames)
		if(any( loc.wildchar >= 0))
		{
			motif_ids <- c();
			for (f in filenames) {
				# no wildchar
				if( loc.wildchar < 0 ) 
				{
					pwm.filenames <- c( pwm.filenames, f );
					next;
				}
				
				wildname <- basename(f);
				fnames = list.files(dirname(f), pattern = glob2rx(wildname), full.names = TRUE, recursive = FALSE)
				
				pwm.filenames <- c( pwm.filenames, fnames );	
				motif_ids <- c( motif_ids, basename( file_path_sans_ext(df$fileNames) ) );
			}
		}
		else
		{
			if( length( motif_ids ) != length( filenames ) )
				stop("The size of 'motif_id' is not equal to 'file.pwm.matrice'.");
			
			pwm.filenames <- filenames;
		}	
		
		for(i in 1:length(filenames)) 
		{
			motif <- try( read.motif( filenames[i], pseudocount= pseudocount, force_even= force_even, ... ) );
			if( class( motif ) !="matrix" )
				stop(paste("Failed to load the pwm file, ", filenames[i], "\n") );

			pwm.matrice[[i]] <- motif;	
		}
		
		tf_info <- data.frame(Motif_ID=motif_ids, TF_Name=NA, Family_Name=NA, DBID=NA);
	}
	else
	{	
		## don't need 'motif_ids' except the format 'pwm.matrix'
		if(!is.null(motif_ids))
			warning("'motif_ids' can not be applied to the formats except pwm.matrix.");

		motif_ids <- NULL;
		format.style <- format;
		
		if(length(format)==1)
		{
			known.format.name <- c("transfac", "jaspar", "mscan", "meme")
			known.format.file <- paste("pwm.format", known.format.name, "txt", sep=".");
		
			if( which(format==known.format.name)>0)
			{
				idx <- which(format==known.format.name);
				known.format.extfile <- system.file("extdata", known.format.file[idx], package="rtfbsdb");

				## Read in the data
				format.style  <- scan(known.format.extfile, what="", sep="\n")
			}
			else
			{
				cat("No default format is decteced, customize format is used.\n")
			}
		}
		
		## Separate elements by one or more whitepace
		format.list <- strsplit(format.style , "[[:space:]]+")
		
		p <- parse_formatted_datafile( format.list, filenames, skip.lines, pseudocount, force_even, ... );
		if(!is.null(p))
		{
			pwm.matrice <- p$pwm;
			tf_info <- p$tf_info;
			pwm.filenames <- p$filename;
		}
		else
		{
			cat("No motifs are imported from the specified file.\n");
			return(tfbs);
		}	
	}
	
	err.missed  <- 0;
	err.empty   <- 0;
	idx.missing <- c();
	idx.new     <- c();

	for ( i in 1:length(pwm.matrice) )
	{
		if( class(pwm.matrice[i])=="try-error" ) 
		{
			# cat("! Can not find PWM file for motif ID=", motif_ids[i], ".\n" );
			err.missed <- err.missed + 1;
			next;
		}
	
		if(NROW(pwm.matrice[i])==0)
		{
			# cat("! No A C G T values in the PWM file for motif ID=", motif_ids[i], ".\n" );
			err.empty <- err.empty + 1;
			next;
		}
		
		## the motif info will be retrieved from 'tf_missing' table
		if( !is.null(tfbs@tf_missing) && ncol(tfbs@tf_missing)!=0 && nrow(tfbs@tf_missing)!=0)
		{
			motif_id <- as.character(tf_info[i, "Motif_ID"]);
			idx0 <- which( motif_id == as.character(tfbs@tf_missing$Motif_ID )) ;
			if( length(idx0) == 1 )
			{
				idx.missing <- c(idx.missing , i); 
				next;
			}
		}

		idx.new <- c(idx.new, i);
	}

	if( err.empty + err.missed > 0 )
		cat("!", err.empty + err.missed, "PWM file(s) are failed to be loaded ( Missing pwm.matrice :", err.missed, ", Empty pwm.matrice :", err.empty, ").\n");

	if(length(idx.missing)>0)
	{
		id.missing           <- as.character( tf_info[idx.missing, "Motif_ID"] );
		nidx                 <- match( id.missing, as.character(tfbs@tf_missing[, "Motif_ID"] ) );
		tf_missing           <- tfbs@tf_missing[ nidx, , drop=F];

		tfbs@tf_info         <- rbind( tfbs@tf_info, tf_missing );
		tfbs@tf_missing      <- tfbs@tf_missing[ -nidx, , drop=F];

		tfbs@pwm             <- append( tfbs@pwm, pwm.matrice[idx.missing] );
		
		tfbs@filename        <- c( tfbs@filename, pwm.filenames[idx.missing] );
		tfbs@mgisymbols      <- c( tfbs@mgisymbols, id.missing ) ; 
		tfbs@ntfs            <- as.integer( length(tfbs@filename) );
		tfbs@distancematrix  <- matrix(, nrow=0, ncol=0);
		tfbs@cluster         <- matrix(, nrow=0, ncol=0); 
		tfbs@expressionlevel <- as.data.frame( NULL );
	}

	if(length(idx.new)>0)
	{
		id.new              <- as.character( tf_info[idx.new, "Motif_ID"] );
		tf_info			    <- tf_info[idx.new, ,drop=F];
		
		## merge two data frame with the same column names as much as possible.
		if (!is.null(tfbs@tf_info) && ncol(tfbs@tf_info)!=0 && nrow(tfbs@tf_info)!=0 )
		{
			rbind.dataframe<- function(df1, df2)
			{
				cols <- unique(c(colnames(df1), colnames(df2)));
				df1 <- lapply( cols, function(x) { if( x %in% colnames(df1)) as.character(df1[, x]) else return(NA);} );
				df2 <- lapply( cols, function(x) { if( x %in% colnames(df2)) as.character(df2[, x]) else return(NA);} );
				
				df1 <- data.frame( do.call( cbind, df1 ), stringsAsFactors=FALSE);
				df2 <- data.frame( do.call( cbind, df2 ), stringsAsFactors=FALSE);

				df <- rbind( df1, df2 );
				colnames(df) <- cols;
				df;
			}

			## merge two (maybe different) data frame 
			tf_info <- rbind.dataframe(tfbs@tf_info, tf_info);
		} 

		tfbs@pwm             <- append( tfbs@pwm, pwm.matrice[idx.new] );
		tfbs@tf_info         <- tf_info;
		tfbs@filename        <- c( tfbs@filename, pwm.filenames[idx.new] );
		tfbs@mgisymbols      <- c( tfbs@mgisymbols, id.new ) ; 
		tfbs@ntfs            <- as.integer( length(tfbs@filename) );
		tfbs@distancematrix  <- matrix(, nrow=0, ncol=0);
		tfbs@cluster         <- matrix(, nrow=0, ncol=0); 
		tfbs@expressionlevel <- as.data.frame( NULL );
	}
	
	cat("*", length(idx.new) + length(idx.missing), " motif(s) have bern imported.\n");
	
	return(tfbs);
}

# format.style including TAGs
# '>' '[' ']' 
# $SKIP n
# $REPEAT
# $LOM        == Line of motif
# $EOM        == End of Motif	    
# $Motif_ID
# $TF_Name
# $A+, $C+, $G+, $T+
# $A $C $G $T
# $Description
# $anyword

## parser one data file with the specific format style
## return all pwm info.

parse_formatted_datafile <- function( format.style, filenames, skip.lines, pseudocount= -7, force_even= FALSE, ... )
{
	merge_pwm_info<-function( pwms )
	{
		tf_info  <- NULL;
		
		pwm.mat <- lapply(pwms, function(pwm) {return(pwm$pwm_mat)});
		pwm.info <- unlist(lapply(pwms, function(pwm) {return(pwm$tf_info)}));
		pwm.motif_id <- unlist(lapply(pwms, function(pwm) {return(pwm$Motif_ID)}));
		pwm.tf_name <- unlist(lapply(pwms, function(pwm) {return(pwm$TF_Name)}));
		pwm.filename <- unlist(lapply(pwms, function(pwm){return(pwm$filename)}));

		if(!is.null(pwm.info))
		{
			tf_tags <- unique( attr(pwm.info, "names") );

			tf_info <- lapply(pwms, function(pwm) { return( pwm$tf_info[tf_tags] )}) 
			tf_info <- do.call(rbind, tf_info);
			colnames( tf_info ) <- tf_tags;
			
			tf_info <- data.frame( Motif_ID = pwm.motif_id, TF_Name = pwm.tf_name, tf_info, stringsAsFactors=FALSE );
		}
		else
			tf_info <- data.frame( Motif_ID = pwm.motif_id, TF_Name = pwm.tf_name, stringsAsFactors=FALSE );
		
		
		return(list(pwm=pwm.mat, tf_info=tf_info, filename=pwm.filename));
	}
	
	make_ACGT_matrix<-function( pwm )
	{
		mat <- NULL;
		
		if( !is.na(pwm$tags["A+"]) )
		{
			VA <- lapply(strsplit(pwm$tags["A+"], ","), as.numeric)
			VC <- lapply(strsplit(pwm$tags["C+"], ","), as.numeric)
			VG <- lapply(strsplit(pwm$tags["G+"], ","), as.numeric)
			VT <- lapply(strsplit(pwm$tags["T+"], ","), as.numeric)
		
			mat <- cbind(A = as.numeric(VA[[1]]), 
						 C = as.numeric(VC[[1]]), 
						 G = as.numeric(VG[[1]]), 
						 T = as.numeric(VT[[1]]) );

		}	
		
		if( !is.na(pwm$tags["A"]) )
		{
		 	VA <- pwm$tags[ which(attr(pwm$tags, "names")=="A")] 
		 	VC <- pwm$tags[ which(attr(pwm$tags, "names")=="C")] 
		 	VG <- pwm$tags[ which(attr(pwm$tags, "names")=="G")] 
		 	VT <- pwm$tags[ which(attr(pwm$tags, "names")=="T")] 

			mat <- cbind(A = as.numeric(VA), 
						 C = as.numeric(VC), 
						 G = as.numeric(VG), 
						 T = as.numeric(VT) );
						 
		}
	
		if(!is.null(mat))
		{
			if( length( which(rowSums(mat)==0) )>0 )
			{
				warning("Invalid row(s) in the PWM matrix is removed");
				mat <- mat[ - which(rowSums(mat)==0), ];
			}
			
			mat <- log(mat/rowSums(mat));
			
			rownames(mat) <- c(1:NROW(mat));

			## Set a pseudocount of ~0.1%.
			mat[ mat == -Inf ] <- pseudocount;

			## Off-by-one bug for odd motifs.  Not sure why?!  For now hack it.
			if((NROW(mat) %% 2) == 1 & force_even) 
				mat <- rbind(mat, log(c(0.25, 0.25, 0.25, 0.25)));
		
		}
		
		return(mat);
	}

	extract_tag_value <- function( pwm, tag_id)
	{
		value <- NULL;
		idx <- which(attr(pwm$tags, "names")==tag_id)
		if(!is.na(idx[1]))
			value <- pwm$tags[ idx[1] ];
			
	 	return( value ); 
	}

	make_tfinfo_dataframe<-function( pwm )
	{
		tag_ids <- attr(pwm$tags, "names");
		
		known_tags <- c("A", "T", "C", "G", "A+", "T+", "C+", "G+", "", "LOM", "Motif_ID", "TF_Name");
		for(tag in known_tags)
			tag_ids <- tag_ids[ tag_ids != tag ];
		
		df <- NULL;
		if ( length(tag_ids)>0 )
		{
			df <- pwm$tags[ match( tag_ids, attr(pwm$tags, "names")) ];
			names(df) <- tag_ids;
		}
		return(df);	
	}	
	
	## remove comments from format style list
	is.comment <- unlist( lapply(format.style, function(s){ substring(s[1],1,1)=="#" }));
	format.style <- format.style[!is.comment];
	
	## split '>TAG' into two parts.
	for(i in 1:length(format.style))
		if( substring(format.style[[i]][1],1,1) == ">" && format.style[[i]][1] != ">" )
			format.style[[i]]<- c( substring(format.style[[i]][1],1,1), substring(format.style[[i]][1],2), format.style[[i]][-1])

	pwms <- list();
	
	for(K in length(filenames))
	{
		lines <- scan( filenames[K], what="", sep="\n" );
		
		## skip some comments 
		if(skip.lines>0)
			lines <- lines[-c(1:skip.lines)];

		lines <- strsplit(lines, "[[:space:]]+")
		i.line <- 1;

		filename <- paste(filenames[K], "#", i.line, sep="")
		cpwm <- list( pwm_mat = NULL, cursor = 1, status = "", next.step = "", filename = filename);  

		while( i.line <= length(lines) )
		{
			cpwm <- parser_aline( cpwm, format.style, lines[[i.line]])
			
			if( cpwm$status=="ERROR" )
			{
				# Something is wrong, skip error and start over.
				cat("\nThe error occured at line ", i.line, "\n" );
				if(i.line-1>0) cat( "   ", i.line-1, "|", lines[[i.line-1]], "\n");
				cat( "***", i.line, "|", lines[[i.line]], "\n");
				if(i.line+1>0) cat( "   ", i.line+1, "|", lines[[i.line+1]], "\n");
				
				# Don't proceed anymore, quit!
				break;
			}

			if( cpwm$status=="EOM" )
			{
				# make PWM matrix @todo;
				cpwm$pwm_mat  <- make_ACGT_matrix( cpwm );
				cpwm$tf_info  <- make_tfinfo_dataframe( cpwm );
				cpwm$Motif_ID <- extract_tag_value( cpwm, "Motif_ID");
				cpwm$TF_Name  <- extract_tag_value( cpwm, "TF_Name");
				cpwm$tags     <- NULL;
#show(cpwm);		
				pwms[[length(pwms)+1]] <- cpwm;
				filename <- paste(filenames[K], "#", i.line, sep="")
				cpwm <- list( pwm_mat = NULL, cursor = 1, status = "", next.step = "NEXT", filename= filename );  
			}

			if( cpwm$next.step == "SKIP" ) 
				i.line <- i.line + cpwm$skip.line
			else if( cpwm$next.step == "NEXT" ) 	
				i.line <- i.line + 1
			else if ( cpwm$next.step == "ROLLBACK")
			{
				#no change the current data line
			}
			else
				cat("UNEXPECTED status=", cpwm$next.step, "\n")
		}
	}

	if(cpwm$status=="ERROR")
		return(NULL)
	else	
		return( merge_pwm_info(pwms));
}

## parser one data line with the current format style
## return the info and status from the data line.

parser_aline<-function( pwm, format.style, data.line)
{
	## put dummy code to pass the check(R CMD check rtfbsdb --as-cran)
	s <- "";

	start_with<-function(s, a ) {return(substring(s,1,1)==a)}
	ending_with<-function(s, a ) {return(substring(s,nchar(s),1)==a)}

	## split the string of '>TAGS' 
	if( substring( data.line[1],1,1 ) == ">" && data.line[1] != ">" )
		data.line<- c( substring( data.line[1],1,1 ), substring( data.line[1],2 ), data.line[-1])

	## data.line: one line from data file
	## bline: split '[TAGS' or 'TAGS]' into two parts

	bline <- c();
	for(i in 1:length(data.line))
	{
		## split the string of '[TAGS' 
		if( start_with( data.line[i], "[") && data.line[i] != "[" )
			bline<- c( bline, substring(data.line[i],1,1), substring(data.line[i],2) )
		## split the string of ']TAGS' 
		else if( start_with( data.line[i], "]") && data.line[i] != "]" )
			bline<- c( bline, substring(data.line[i],1,1), substring(data.line[i],2) )
		## split the string of 'TAGS]' 
		else if( ending_with( data.line[i], "]")  && data.line[i] != "]")
			bline<- c( bline, substring(data.line[i], 1, nchar(data.line[i])-1 ), substring(data.line[i], nchar(s), 1) )
		else if (nchar(data.line[i])>0)
			bline <- c( bline, data.line[i] );
	}

#cat("Data=", data.line	, "\n");
	r.match <- list( matched = FALSE );
	
	if( pwm$cursor <= length(format.style) )
	{
		r.match <- match_aline( format.style[[pwm$cursor]], bline);
#cat("Rule=", format.style[[pwm$cursor]], "\n");
		if( !r.match$matched ) 
			## if the current format does not match the data, move cursor to the next line of format data.
			if( pwm$cursor < length(format.style) )
			{
				r.match <- match_aline( format.style[[pwm$cursor+1]], bline);
#cat("Rule=", format.style[[pwm$cursor+1]], "\n");
				if(r.match$matched) 
					pwm$cursor <- pwm$cursor + 1;
			}
	}
	else
	{
		pwm$status <- "ERROR";
		return(pwm);
	}	
	
	if( r.match$matched )
	{
		if( length(r.match$values) > 0 )
		{
			## put the data into pwm	
			tags <- r.match$values; 
			## remove dollar symbols
			names(tags) <- substring( r.match$names, 2 ); 
			## merge tags 
			pwm$tags <- c(pwm$tags, tags);
		}

		pwm$status    <- r.match$status;
		pwm$next.step <- "NEXT";
		
		## '$REPEAT'
		pwm$is.repeat <- r.match$repeat.word;
		
		## '$SKIP n' (n>1)
		if ( r.match$skip.line > 1)
		{
			pwm$next.step <- "SKIP";
			pwm$skip.line <- r.match$skip.line - 1; 
		}
		
		## move cursor if not '$REPEAT'
		if( !r.match$repeat.keyword )	
		{
			pwm$cursor <- pwm$cursor + 1;
			if( pwm$cursor > length(format.style) )
				pwm$status <- "EOM";
		}		
	}
	else
	{
		if( pwm$cursor >= length(format.style)  )
		{	
			## try to match the current line with the first format data.
			r.first <- match_aline( format.style[[1]], bline);

			## if successful, it's a new data line
			if( r.first$matched )
			{
				pwm$status <- "EOM";
				pwm$next.step <- "ROLLBACK";
			}
		}
		else
			pwm$status <- "ERROR";
	}

	return(pwm);
}


## parser one data line with the current format line
## return the info from the data line.

match_aline <- function( format.line, data.line)
{	
	## put dummy code to pass the check(R CMD check rtfbsdb --as-cran)
	s <- "";
	
	start_dollar<-function(s) { return(substring(s,1,1)=="$")}
	is_numeric <- function(s) { suppressWarnings(!is.na(as.numeric(s))) };
	
	## if no tag at the head of line, the first word must atch the format data.
	if( !start_dollar( format.line[1] ) && data.line[1] != format.line[1] )
		return(list(matched=FALSE));
	
	## f: the word cursor of format line
	## b: the word cursor of data line
	
	f <- 1; b <- 1;
	skip.line <- 0;
	repeat.keyword <- FALSE;
	status <- "TAG"; 
	matched <- FALSE;

	tag.values <- c();
	tag.names <- c();
	tag.plus <- c();

	while( b <= length(data.line))
	{	
		## word macthed 
		if( data.line[b]==format.line[f] && 
			!start_dollar(format.line[f]) &&
			!start_dollar(data.line[b]) )
		{
			b <- b + 1;
			f <- f + 1;
			matched = TRUE;
			next;
		}

		## '$SKIP n'
		if ( format.line[f]=="$SKIP" )
		{
			status <- "SKIP";
			skip.line <- 1;
			matched <- TRUE;
			if(!is.na(format.line[f+1]) && is_numeric( format.line[f+1] ) )
			{
				## Here 'skip.line' just original value
				skip.line <- as.numeric(format.line[f+1]);
			}	
			
			break;
		}

		## '$REPEAT'
		if ( format.line[f]=="$REPEAT" )
		{	
			repeat.keyword <- TRUE;
			break;
		}	

		## '$LOM'
		if ( format.line[f]=="$LOM"  )
		{
			if( is_numeric(data.line[b]) )
			{
				tag.values <- c(tag.values, data.line[b]);
				tag.names <- c(tag.names, format.line[f]);
				
				b <- b + 1;
				f <- f + 1;
				status <- "ACGT"; 
				matched <- TRUE;
				next;
			}
			else
			{
				matched <- FALSE;
				break;
			}	
		}

		## '$EOM'
		if ( format.line[f]=="$EOM" )
		{
			if(is.na(data.line[b]))
			{
				status <- "EOM"; 
				matched <- TRUE;
				break;
			}
			else
			{
				matched <- FALSE;
				break;
			}
		}

		## '$ACGT'
		if ( format.line[f] %in% c("$A", "$C", "$G", "$T" )  )
		{
			if( is_numeric(data.line[b]) )
			{
				tag.values <- c(tag.values, data.line[b]);
				tag.names <- c(tag.names, format.line[f]);

				b <- b + 1;
				f <- f + 1;
				status <- "ACGT"; 
				matched <- TRUE;

				next;
			}
			else
			{
				matched <- FALSE;
				break;
			}
		}

		## '$ACGT+'
		if ( format.line[f] %in% c("$A+", "$C+", "$G+", "$T+" ) )
		{
			if(is_numeric(data.line[b]))
			{
				if(length(tag.plus)==0) tag.plus <- format.line[f] ;
				tag.plus <- c(tag.plus, data.line[b]);

				b <- b + 1;
				status <- "ACGT+"; 
				matched <- TRUE;
				next;
			}
			else
			{
				f <- f + 1;
				break;
			}
		}
		
		## other $TAGs
		if ( start_dollar(format.line[f]))
		{
			tag.values <- c(tag.values, data.line[b]);
			tag.names <- c(tag.names, format.line[f]);
			
			b <- b + 1;
			f <- f + 1;
			matched <- TRUE;
			status <- "TAG";
			next;
		}
	}
	
	if ( b > length(data.line) && f <= length(format.line) && matched )
	{
		if ( format.line[f]=="$EOM" )
		{
			status <- "EOM"; 
			matched <- TRUE;
		}

		if (format.line[f]=="$SKIP" )
		{
			status <- "SKIP";
			if(!is.na(format.line[f+1]) && is_numeric( format.line[f+1] ))
				## Here 'skip.line' just original value
				skip.line <- as.numeric(format.line[f+1]);
		}

		if (format.line[f]=="$REPEAT" )
		{	
			repeat.keyword <- TRUE;
		}	
	}
	
	if( length( tag.plus)>0)
	{
		tag.values <- c( tag.values, paste( tag.plus[-1], collapse=",", sep=" ") );
		tag.names  <- tag.plus[1];
	}
	
	return(list(matched=matched, status=status, values=tag.values, names=tag.names, repeat.keyword=repeat.keyword, skip.line=skip.line))
}
