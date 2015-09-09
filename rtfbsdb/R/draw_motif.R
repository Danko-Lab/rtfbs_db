#' Basic functions for drawing TFs.
#'
#' tfbs.drawLogo
#'
#' @file.pdf
#' @tf_id
#' @motif_id
#' @tf_name
#' @family_name
#' @tf_status
#' @group_by

tfbs_drawLogo <- function(tfbs, file.pdf=NULL, index=NULL, tf_id=NULL, motif_id=NULL, tf_name=NULL, family_name=NULL, tf_status=NULL, groupby=NULL) 
{
	idx.select <- c();

	if( !is.null(index) ) idx.select <- c(idx.select, index);
	if( NROW( tfbs@tf_info )>0  && !is.null(tf_id) ) 
		idx.select <- c( idx.select, which(!is.na(match(tfbs@tf_info$TF_ID, tf_id ))) );
	if( NROW( tfbs@tf_info )>0  && !is.null(motif_id) ) 
		idx.select <- c( idx.select, which(!is.na(match(tfbs@tf_info$Motif_ID, motif_id ))) );
	if( NROW( tfbs@tf_info )>0  && !is.null(tf_name) ) 
		idx.select <- c( idx.select, which(!is.na(match(tfbs@tf_info$TF_Name, tf_name ))) );
	if( NROW( tfbs@tf_info )>0  && !is.null(family_name) ) 
		idx.select <- c( idx.select, which(!is.na(match(tfbs@tf_info$Family_Name, family_name ))) );
	if( NROW( tfbs@tf_info )>0  && !is.null(tf_status) ) 
		idx.select <- c( idx.select, which(!is.na(match(tfbs@tf_info$TF_Status, tf_status ))) );
	if( is.null(index) && is.null(tf_id) && is.null(motif_id) && is.null(tf_name) && is.null(family_name) && is.null(tf_status)) 
		idx.select <- c(1:tfbs@ntfs);

	if( length(which(is.na(idx.select)))>0 )
	{
		cat("!", length(which(is.na(idx.select))), "motifs can not be found in the tfbs object.");
		idx.select <- idx.select[ !is.na(idx.select) ];
	}

	idx.select <- sort(unique(idx.select));

	draw_viewport <- function(i, xaxis = TRUE, yaxis = TRUE, cex = 1 )
	{
		tf_name <- tfbs@mgisymbols[i];	
		if ( NROW(tfbs@tf_info) > 0 )
			tf_name <- paste( tfbs@tf_info[ i, "TF_Name"], " (Motif_ID:", tfbs@tf_info[ i, "Motif_ID"], "/DBID:", tfbs@tf_info[ i, "DBID"], ")", sep="");

		pushViewport( viewport(x=0, y=0.94, width=1, height=0.05, just=c("left","bottom")) );
		grid.text( tf_name, rot=0, gp=gpar(cex=cex), check.overlap=T);  
		popViewport();

		pushViewport( viewport(x=0, y=0, width=1, height=0.94, just=c("left","bottom")) );
		seqLogo( exp(t(tfbs@pwm[[i]])), xaxis = xaxis, yaxis = yaxis)
		popViewport();
	}
	
	tfbs.grpupdby<-c("Family_Name", "TF_Name", "TF_Status", "Motif_Type");
	if( !is.null(groupby) && length( which(groupby == tfbs.grpupdby ) ) == 0)
	{
		cat("The availabe value for groupby are Family_Name, TF_Name, TF_Status and Motif_Type.\n");
		groupby <- NULL;
	}	 

	if( !is.null(file.pdf) && !is.na(file.pdf) ) pdf(file.pdf); 

	if(is.null(groupby))	
	{
		for(i.motif in idx.select)
		{
			if( i.motif<=0 || i.motif>tfbs@ntfs)
				next;

			grid.newpage()
			vp1 <- viewport(x=0, y=0, width=1, height=1, just=c("left","bottom"))
			pushViewport(vp1);
			draw_viewport(i.motif);
			popViewport();
		}
	}
	else
	{
		groups <- unique( as.character( tfbs@tf_info[idx.select, c(groupby)] ) )

		for(k in 1:length(groups) )
		{
			idx.page <- idx.select[ which( tfbs@tf_info[idx.select, c(groupby)] == groups[k])];

			grid.newpage();
			vp1 <- viewport(x=0, y=0, width=1, height=1, just=c("left","bottom"))
			pushViewport(vp1);

			for(i in 1:length(idx.page) )
			{
				if(i>10 && i%%10==1)
				{
					popViewport();
					grid.newpage();
					vp1 <- viewport(x=0, y=0, width=1, height=1, just=c("left","bottom"))
					pushViewport(vp1);
				}
				
				i.motif <- idx.page[i];
				if( i.motif<=0 || i.motif>tfbs@ntfs )
					next;
		
				vp1 <- viewport(y = 0.8 - (((i-1)%%10)%/%2)/5, x=(i+1)%%2/2, width=0.5, height=0.2, just=c("left","bottom"))
				pushViewport(vp1);
				draw_viewport(i.motif, xaxis = FALSE, yaxis = FALSE, cex=0.6 );
				popViewport();

			}

			popViewport();
		}
	}

	if( !is.null(file.pdf) && !is.na(file.pdf) ) dev.off(); 
}

#' drawing TFs for each clustering.
#'
#' tfbs.drawLogosForClusters
#'
#' @file.pdf
#' @cluster.mat

tfbs_drawLogosForClusters <- function(tfbs,  file.pdf=NULL) 
{
	if( NROW( tfbs@cluster )  == 0)
		stop("No cluster information in the tfbs object.")
	
	if( !is.null(file.pdf) && !is.na(file.pdf) )
		if( !check_folder_writable( file.pdf ) ) 
			stop( paste("Can not create pdf file: ",file.pdf ) );

	if( !is.null(file.pdf) && !is.na(file.pdf) )
	{
		r.try <- try ( pdf( file.pdf ) );	
		if(class(r.try)=="try-error")
			stop("Failed to create PDF file for motif logos.\n");
	}
		
	draw_tf_name <- function(tfbs, idx, nmax.motifs)
	{
		tf_name <- tfbs@mgisymbols[idx];	
		if ( NROW(tfbs@tf_info) > 0 )
			tf_name <- tfbs@tf_info[ idx, "TF_Name"];
			
		n.block  <- round(nmax.motifs/2);
		pushViewport( viewport(x=0.01, y=0, width=0.05, height=1, just=c("left","bottom")) );
		  	
		str.inches <- strwidth(tf_name, units = 'in');
		y.scale <- convertHeight(unit(str.inches, "inches"), "native", valueOnly=T );
		
		cex <- 1;
		if( y.scale> 0.8 ) cex <- 0.8/y.scale;
		if( cex > 1 ) cex <- 1;

		grid.text( tf_name, rot=90, gp=gpar(cex=cex), check.overlap=T);  

		popViewport();
	}

	cluster.mat <- tfbs@cluster;

	for(i in 1:max(cluster.mat[,2])) 
	{
		motifs_in_cluster <- cluster.mat[which(cluster.mat[,2]==i), 1];
		nmotifs <- length(motifs_in_cluster);
			
		nmax.motifs <- nmotifs;
		# at least4 motif space will be designed in each page.
		if( nmax.motifs <= 4 ) nmax.motifs <- 4;
			
		grid.newpage();
		for(j in seq(1, nmax.motifs, 2)) {
			
			vp1 <- viewport(x=0, y=((j-1)/2)/ceiling(nmax.motifs/2), width=0.5, height=1/ceiling(nmax.motifs/2), just=c("left","bottom"));
			pushViewport(vp1); 

			#seqLogo(makePWM(exp(t(tfbs@pwm[[motifs_in_cluster[j]]]))), xaxis = FALSE, yaxis = FALSE);
			if ( j <= nmotifs )
			{
				draw_tf_name(tfbs, motifs_in_cluster[j], nmax.motifs);
				pushViewport( viewport(x=0.04, y=0, width=0.96, height=1, just=c("left","bottom")) );
				seqLogo(exp(t(tfbs@pwm[[motifs_in_cluster[j]]])), xaxis = FALSE, yaxis = FALSE);
				popViewport();
			}
			popViewport();
			  
			vp1 <- viewport(x=0.5, y=((j-1)/2)/ceiling(nmax.motifs/2), width=0.5, height=1/ceiling(nmax.motifs/2), just=c("left","bottom"));
			pushViewport(vp1); 
			#seqLogo(makePWM(exp(t(tfbs@pwm[[motifs_in_cluster[j]]]))), xaxis = FALSE, yaxis = FALSE);
			if ( j+1 <= nmotifs )
			{
				draw_tf_name(tfbs, motifs_in_cluster[j+1], nmax.motifs);
				pushViewport( viewport(x=0.04, y=0, width=0.96, height=1, just=c("left","bottom")) );
				seqLogo( exp(t(tfbs@pwm[[motifs_in_cluster[j+1]]])), xaxis = FALSE, yaxis = FALSE);
				popViewport();
			}
			
			popViewport();
		}
	}

	if( !is.null(file.pdf) && !is.na(file.pdf) )
		dev.off();	
}