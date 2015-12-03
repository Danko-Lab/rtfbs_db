output_motif_report<-function( tfbs, df.motif, file.pdf, report.size, report.title, report.style, report.note="")
{
	scheme1 <- data.frame( min.col = c( 0,0,1 ), max.col = c( 0.9, 0.9, 0.9 ) );
	scheme2 <- data.frame( min.col = c( 0,1,0 ), max.col = c( 0.6, 0.3, 0.6 ) );
	
	get_rgbcol<-function(pvalue, p.min=0.001, p.max=1, log10=T, scheme=1)
	{
		if(pvalue > p.max) pvalue <- p.max;
		if(pvalue < p.min) pvalue <- p.min;
		
		if(log10)
		{
			pvalue <- log10(pvalue)
			p.min  <- log10(p.min)
			p.max  <- log10(p.max)
		}
		
		if(scheme==1) col.scheme <- scheme1 else col.scheme <- scheme2;
		
		c1 <- (col.scheme$min.col[1] - col.scheme$max.col[1]) / (p.min-p.max)*(pvalue-p.max) + col.scheme$max.col[1];
		c2 <- (col.scheme$min.col[2] - col.scheme$max.col[2]) / (p.min-p.max)*(pvalue-p.max) + col.scheme$max.col[2];
		c3 <- (col.scheme$min.col[3] - col.scheme$max.col[3]) / (p.min-p.max)*(pvalue-p.max) + col.scheme$max.col[3];

		c <- rgb(c1,c2,c3)

		return(c);
	}

	drawlegend<-function( x0, y0, width, height, title, pv.min=1e-6, pv.max=1, pv.log10=T, scheme=1)
	{
		bar.len <- 50;
		if(pv.log10)
			pv <- 10^seq(log10(pv.min), log10(pv.max), length.out=bar.len)
		else
			pv <- seq(pv.min, pv.max, length.out=bar.len); 
		
		bar.width <- width*0.6;
		for( i in 1:length(pv) )
		{
			col <- get_rgbcol( pv[i], pv.min, pv.max, log10=pv.log10, scheme=scheme);
			grid.rect(  x0 + 0.3*width + i*bar.width/bar.len, 
						y0, 
						width = bar.width/bar.len, 
						height = height, 
						gp = gpar(fill=col, col=NA), 
						just = c("right", "centre") );
		}
		
		if( pv.log10 )
		{
			p1<- (log10(pv.max) - log10(0.05) )*bar.width/(log10(pv.max) - log10(pv.min) );
	
			grid.text(  ".05", 
						p1 + x0 + 0.3*width, 
						y0 - height, 
						gp = gpar(cex=0.4), 
						just = c("centre",  "centre"));
		}
		else
		{
			p1 <- ( 1 - pv.min) * bar.width/(pv.max - pv.min );	

			grid.text(1, p1+x0+0.3*width,  y0-height, gp=gpar(cex=0.4), just=c("centre",  "centre"));
		}
		
		grid.lines( x = c( p1 + x0 + 0.3*width, p1 + x0 + 0.3*width), 
					y = c( y0 - height/2, y0 + height/2), 
					gp = gpar(col="white", fill = "white", lwd = 0.01) );
		
		grid.text(  title,  
					x0 + 0.25*width,  
					y0, 
					gp = gpar(cex=0.4), 
					just = c("right",  "centre"));
		grid.text( sprintf("%1g", pv.min), 
					x0 + 0.28*width,  
					y0 - height, 
					gp = gpar(cex=0.4), 
					just = c("centre", "centre"));
		grid.text( sprintf("%1g", pv.max), 
					x0 + 0.9*width,  
					y0 - height, 
					gp = gpar(cex=0.4), 
					just = c("centre", "centre"));
	}

	get_short_value<-function(val, log=T)
	{
		if(log)
		{
			n.log <- ceiling(log10(val));
			if(n.log > -2)
				str <- round(val, digits=2)
			else
				str <- sprintf("<%.1g", val);
			if(val>=0.05) str<-"";
			return(str);
		}
		else
			return(round(val, digits=1));
	}
	
	draw_top <- function()
	{
		grid.text( report.title, x=0.5, y= 0.6,  gp=gpar(cex=0.8), just=c("centre", "centre"));

		grid.lines( x = c(0,1), y=c(0.4,0.4), gp=gpar(col="black", lwd=1) );

		for( k in 1:NROW(report.style) )
		{
			if(as.character(report.style$hjust[k])=="left")
				grid.text(  report.style$header[k], 
							x = report.style$position[k]+0, 
							y = 0.2, 
							rot = 0, 
							gp = gpar(cex=0.5), 
							just=c("left", "centre") ); 
			if(as.character(report.style$hjust[k])=="centre")
				grid.text(  report.style$header[k], 
							x = report.style$position[k]+report.style$width[k]/2, 
							y = 0.2, 
							rot = 0, 
							gp = gpar(cex=0.5), 
							just = c("centre", "centre") ); 
			if(as.character(report.style$hjust[k])=="right")
				grid.text( report.style$header[k], 
							x = report.style$position[k]+report.style$width[k], 
							y = 0.2, 
							rot = 0, 
							gp = gpar(cex=0.5),  
							just = c("right", "centre") ); 
		}
		
		grid.lines(x = c(0,1), y=c(0, 0), gp=gpar(col="black", lwd=1) );
	}

	draw_bottom <- function()
	{
		grid.text(report.note, x=1, y=0.95,  gp=gpar(cex=0.4), just=c("right", "top"));

		x0 <- 0.8;
		for( k in 1:NROW( report.style) )
		{
			if(report.style$style[k]=="bar")
			{
				p.min  = as.numeric(as.character(report.style$extra1[k])); 
				p.max  = as.numeric(as.character(report.style$extra2[k])); 
				log10  = as.numeric(as.character(report.style$extra3[k]));
				scheme = as.numeric(as.character(report.style$extra4[k]));
					
				drawlegend( x0, 0.5, width=0.2, height=0.3, title=report.style$header[k], p.min, p.max, log10, scheme );
				x0 <- x0 - 0.2;
			}										
		}
	}

	draw_item<-function(r.comp.sort, idx.start, idx.stop)
	{
		M40 <- view.lines;
	
		y0 <- 0;
		for(i in idx.start:idx.stop)
		{
			y0 <- 1 - (i - idx.start +1 )/M40;
			for(k in 1:NCOL(df.motif))
			{
				if(report.style$style[k]=="text")
				{
					if(as.character(report.style$hjust[k]) == "left")
						grid.text(  df.motif[i, k], 
									x = report.style$position[k]+0, 
									y = y0+0.005, 
									rot = 0, 
									gp = gpar(cex=0.5), 
									check.overlap = F, 
									just = c("left", "centre") ); 
					if(as.character(report.style$hjust[k]) == "centre")
						grid.text( df.motif[i, k], 
									x = report.style$position[k]+report.style$width[k]/2, 
									y = y0+0.005, 
									rot = 0, 
									gp = gpar(cex=0.5), 
									check.overlap = F, 
									just=c("centre", "centre") ); 
					if(as.character(report.style$hjust[k])=="right")
						grid.text( df.motif[i, k], 
									x = report.style$position[k]+report.style$width[k], 
									y = y0+0.005, 
									rot = 0, 
									gp = gpar(cex=0.5), 
									check.overlap = F, 
									just = c("right", "centre") ); 
				}						

				if(report.style$style[k]=="bar")
				{			
					col.fill <- get_rgbcol( df.motif[i,k], 
										p.min = as.numeric(as.character(report.style$extra1[k])), 
										p.max = as.numeric(as.character(report.style$extra2[k])), 
										log10 = as.numeric(as.character(report.style$extra3[k])),
										scheme = as.numeric(as.character(report.style$extra4[k]))) ; 
					grid.rect(  x = report.style$position[k] + report.style$width[k]/2.0, 
								y = y0+0.005, 
								width = report.style$width[k]-0.02, 
								height = 1/M40-0.006, 
								gp = gpar(lwd=0.01, 
								fill = col.fill) );

					str.val <- get_short_value( df.motif[i,k], as.numeric(as.character(report.style$extra3[k])) );

					grid.text(  str.val, 
								x = report.style$position[k] + report.style$width[k]/2.0, 
								y = y0+0.005, 
								rot = 0, 
								gp = gpar( cex=0.4, col="white"), 
								check.overlap = F ); 
				}
				
				if(report.style$style[k]=="logo" && !is.null(tfbs@tf_info))
				{			
					pushViewport( viewport( y = y0, 
								x = report.style$position[k], 
								width = report.style$width[k], 
								height = 1/M40, 
								just = c("left","bottom")));

					idx <-which( as.character(df.motif[i,k]) == as.character(tfbs@tf_info$Motif_ID) );
					if(length(idx)>1)
					{
						warning(paste("Multiple matrice in the tfbs object for Motif ID:", as.character(df.motif[i,k]), 
							", first matrix (index:", idx[1], ") is used to draw logo.\n"));
						idx <- idx[1];
					}
					
					seqLogo( exp(t(tfbs@pwm[[idx]])), xaxis = FALSE, yaxis = FALSE);
					popViewport();
				}
				
			}
		}
	
		grid.lines( x = c(0, 1), y=c( y0-1/M40, y0-1/M40 ), gp=gpar(col="black", lwd=1) );
		return(1 - ( y0 - 1/M40) );
	}

	view.height <- 1;
	view.lines  <- 28;
	if(!is.na(file.pdf)) 
	{
		if (report.size=="letter") 
		{
			view.height <- 10/6;
			view.lines  <- 28*10/6;
			pdf( file.pdf, width=8, height=10 );
		}
		else
			pdf( file.pdf, width=6, height=6 );

	}

	for(i in 1:ceiling(NROW(df.motif)/view.lines))  
	{
		grid.newpage();
		pushViewport( viewport(x=0, y=0, 
						width=1, height=1, 
						just=c("left","bottom")) );
		
		pushViewport( viewport(x=0.15, y=0.97, 
						width=0.7, height=0.07, 
						just=c("left","top"), 
						xscale = c(0, 1), 
						yscale = c(0, 1) ) )
		draw_top();
		popViewport();

		pushViewport( viewport(x=0.15, y=0.90, 
						width=0.7, height=0.82, 
						just=c("left","top"), 
						xscale = c(0, 1), 
						yscale = c(0, 1) ) )
		y0 <- draw_item( df.motif, (i-1)*view.lines+1, min( i*view.lines, NROW(df.motif) ) );
		popViewport();

		pushViewport( viewport(x=0.15, y=0.9-0.82*y0, 
						width=0.7, height=0.04, 
						just=c("left","top"), 
						xscale = c(0, 1), 
						yscale = c(0, 1) ) )
		draw_bottom();
		popViewport();

	}    

	popViewport();
	
	if(!is.na(file.pdf)) dev.off(); 
}