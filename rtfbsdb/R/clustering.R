tfbs_clusterMotifs <- function(tfbs, subset=NA, pdf.heatmap=NA, method=c("agne", "cors"), group.k=NA)
{
	mat <- tfbs@distancematrix;
	if (is.null(mat) || all(mat==0))
		stop("The distance matrix is not available, run tfbs.getDistanceMatrix firstly.");
	
	if(!missing(subset)) 
		mat <- tfbs@distancematrix[subset,subset,drop=F]
	else
		subset <- 1:NROW(mat);
	
	if(missing(method)) 
		method <- "agne";

	if( !is.na( pdf.heatmap ) )
    	if( !check_folder_writable( pdf.heatmap ) ) 
  		    cat("! Can not create pdf file: ", pdf.heatmap, "\n");
		
	if(method=="cors")
	{
		opt_list <- tfbs_corclustering_bic_optim(mat);

		bics <- lapply(opt_list, function(l)return(l$BIC) );
		min.grps <- which.min( bics);
		clusters <- opt_list[[min.grps]]$clusters;

cat("MIN CLUSRERS=", min.grps, max( opt_list[[min.grps]]$clusters))	
		
		if(!is.na( pdf.heatmap ))
		{
			r.try <- try ( pdf( pdf.heatmap ) );
			if( class(r.try) == "try-error")
				cat("! Failed to write PDF file:", pdf.heatmap, "\n")
			else
			{
				tfbs_drawheatmapForClusters(tfbs, mat, clusters);
				dev.off(); 
			}
		}
	}
	else
	{
		if(missing(group.k))
			group.k <- round( NROW(mat)/3 );
			
		hc1 <- agnes(as.dist((1-mat)^5), diss=TRUE)
		clusters <- cutree(hc1, k=group.k)

		if(!is.na( pdf.heatmap ))
		{
			r.try <- try ( pdf( pdf.heatmap ) );
			if( class(r.try) == "try-error")
				cat("! Failed to write PDF file:", pdf.heatmap, "\n")
			else
			{
           		hc1 <- as.dendrogram(hc1)
           		ord.hc1 <- order.dendrogram(hc1)
           		hc2 <- reorder(hc1, mat[ord.hc1])
           		ord.hc2 <- order.dendrogram(hc2)

           		pal100 <- c("#7E291B","#66E52C","#8F66F0","#58DBE8","#396526","#EAABC1","#E1C33C","#3E3668","#EB3F90","#C3E6A8","#E74618","#66A2E9","#3E7774","#DF9056","#3C2C21","#DF40D7","#6CEF92","#8C5A6B","#BC8AE2","#A03B99","#56AC2D","#389C6C","#E26B7E","#706B4B","#D2E374","#A0A560","#7B1C3E","#49F7DB","#C8C6E9","#414FA2","#95A590","#8A669A","#98A62F","#9E792C","#D69489","#547FE5","#DF6340","#849BAB","#E63A61","#8A386D","#DAE338","#263715","#BBDEE3","#3F1324","#A1E03B","#383544","#76C2E8","#794D38","#DD74E2","#D7BF8E","#E366B9","#894D19","#D57221","#D9B764","#B0303D","#D6BFBC","#757A28","#D991CC","#356344","#E3A22D","#223C36","#83B664","#D8E4C5","#AE9ED9","#37576D","#CF7398","#6BEAB2","#6C9266","#B33662","#4B340E","#57E65B","#E23635","#9B464B","#757089","#578BBA","#A6311B","#B2714E","#457F20","#4EA3B1","#B93286","#73D463","#531914","#6B78C2","#E07467","#8D786E","#515018","#361E40","#AA4FCC","#90D2C4","#71469B","#419C4C","#37558A","#B393B0","#9BE090","#856BDA","#66B593","#47CA84","#5D205B","#672F3F","#59D8C2")
           		pal500 <- rep(pal100, 5)

		   		cuth <- clusters[ ord.hc2 ];

		   		fill.col<- rep("white", length(ord.hc2));
		   		fill.col[ord.hc2] <- pal500[ clusters ];

           		#pl <- levelplot((mat)[ord.hc2, ord.hc2], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="",
           		print( levelplot((mat)[ord.hc2, ord.hc2], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="",
           			colorkey = list(space="left", labels=list(cex=1.5)), 
           			legend = list(
           			    top = list(fun = dendrogramGrob,
           			    args = list(x = hc2, ord = ord.hc2, side = "top", #lwd=2,
           			    size = 7, size.add = 0.5, 
           			    add = list(rect = list(col = "transparent", fill = fill.col )),
           			    type = "rectangle")))) );

		   		dev.off();               
		   }
      }
	}
	
	return( cbind(subset, clusters) );
}


tfbs_corclustering_bic_optim<-function( obs_mat )
{
	get_LR2<-function(obs_mat, clustering_mat, clusters)
	{
		N.sample <- NROW(obs_mat);
		
		v1 <- rep( clusters, NROW(obs_mat));
		v2 <- rep( clusters, each=NROW(obs_mat));

		idxs <- which(v1==v2);

		dif_vec <- c(clustering_mat)[idxs] -  c(obs_mat)[idxs] ;
		sd_mat  <- sd (dif_vec);
		mu_mat  <- mean( dif_vec );

		LR <- sum( d_cor_dist_log( c( clustering_mat )[idxs], c(obs_mat)[idxs], N.sample))  
		
		return(LR);
	}

	get_LR<-function(obs_mat, clustering_mat, clusters)
	{
		N.sample <- NROW(obs_mat);
		
		dif_vec <- c(clustering_mat) -  c(obs_mat);
		sd_mat  <- sd (dif_vec);
		mu_mat  <- mean( dif_vec );

		LR2 <- sum( d_cor_dist_log( c( clustering_mat ), c(obs_mat), N.sample ) ) - N.sample * d_cor_dist_log( 1, 1, N.sample );
		return(LR2/2);
	}
		

	set.peak.range <- 1;
	set.var.range  <- c(  1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90 )

	optim_list <- c();

#show(obs_mat);

	loop <-1;		
	for( vr in set.var.range )
	{
		r_cluster <- list( label=loop, var.range=vr, peak.range = set.peak.range, correlation=NA, clusters=NA, AIC=NA, BIC=NA)
		r_cluster$clusters <- tfbs_corclustering_core( obs_mat, var.range=vr, peak.range = set.peak.range )

		r_cluster$cluster_cor_mat <- clusterCorr_fast( obs_mat, r_cluster$clusters );
			
		n.grp <- max(r_cluster$clusters);

		LR <- get_LR( obs_mat, r_cluster$cluster_cor_mat, r_cluster$clusters );
		r_cluster$AIC <- -2*LR + n.grp*(n.grp-1)/2 * 2 ;
		r_cluster$BIC <- -2*LR + n.grp*(n.grp-1)/2 * log( NROW(obs_mat)*(NROW(obs_mat)-1)/2 ) ;
			
		optim_list[[loop]] <- r_cluster;

		cat("loop=", loop, "/", vr, max(r_cluster$clusters), "LR=",	LR, r_cluster$AIC,  r_cluster$BIC, "\n");

		loop <- loop+1;		
	}
	
	if(0)
	{
		pdf("test-AICs.pdf");
		AICs <- lapply(optim_list, function(l)return(l$AIC) );
		plot( unlist(AICs), type = "b" )
		dev.off();

		pdf("test-BICs.pdf");
		BICs <- lapply(optim_list, function(l)return(l$BIC) );
		plot( unlist(BICs), type = "b" )
		dev.off();
	}
	
	return(optim_list)
}

tfbs_corclustering_core<-function( obs_mat, var.range=0.95, peak.range=0.9 )
{
	for(i in 1:NROW(obs_mat)) obs_mat[i,i] <- 0;

	mat <- obs_mat+1;

	find.peaks <- function( mat.row )
	{
		return( length(which(mat.row >= peak.val ) ) );
	}

	find.correlate.close<- function()
	{
	  pool.idx <- c(1:NROW(mat));	
	  rem.idx <- c();

	  for(i in 1:NROW(mat) )
		 if( max(mat[i,]) < max(mat) * var.range ) rem.idx <- c(rem.idx, i);
		 
	  mat.rem <- mat;	 
	  if(length( rem.idx)>0) 	 
	  {
	  	mat.rem <- mat.rem[-rem.idx, -rem.idx,drop = F];
	    pool.idx <- pool.idx[-rem.idx];
	  }
	  
	  sel.idx <- which(mat.rem == max(mat.rem), arr.ind = TRUE)[1,];

	  loop.add <- 1
	  while(loop.add>0)
	  {
	  	  loop.add <- 0;
		  for(i in 1:NROW(mat.rem) )
		  {
		     if( i %in% sel.idx) next;

			 if(max(mat.rem)>1)
	     	 	if ( (max(mat.rem[i,]) -1 )/( max(mat.rem) - 1 ) < var.range ) next;
	     	
			 mat0 <- mat.rem[ c(sel.idx, i), c(sel.idx, i) ];

			 for(k in 1:NROW(mat0)) mat0[k,k] <- max(mat0);
			 
			 if( min(mat0) != max(mat.rem) )
				 if( (min(mat0)-1)/(max(mat.rem)-1) < var.range ) next;

			 r.ix.unsel <- mat.rem[ i, -sel.idx ] ;
			 r.ix.sel <- mat.rem[ i, sel.idx ] ;


			 sel.idx <- c(sel.idx, i );
	  	     loop.add <- loop.add + 1;
		  }
	  }
	  
	  sel.idx <- pool.idx[ sel.idx ];

	  return(unique(sel.idx));
	}

	grp<- list()
	subset.id <- 1:NROW(mat);
	peak.val <- 1;
	loop <- 1;

	while(NROW(mat)>1)
	{
	   	peak.val  <- max( mat - 1 ) * peak.range + 1;
	   	peak.num  <- apply( mat, 1, find.peaks );
	   	
		inner.idx <- find.correlate.close();
	
	   	grp[[ loop ]]  <- subset.id [ inner.idx ];
	   	subset.id <- subset.id[ -inner.idx ];

		mat0 <- mat[ inner.idx, inner.idx, drop=F ];
	   	mat <- mat[ -inner.idx, -inner.idx ];

		for( mm in 1:NROW(mat0))  mat0[mm,mm]<-max(mat0);

	   	loop <- loop + 1;
	}

	if(NROW(mat)==1)
	{
		grp[[loop]] <- c ( subset.id );
	}

	clusters <- rep( 0, NROW(obs_mat) );
	for(i in 1:length( grp ) )
	{
		clusters[c(grp[[i]])] <- i;
	}

	return(clusters);
}

d_cor_dist_F<-function(a, b, c, z, N)
{
	r  <- 1
	rn <- a * b * z / c / 1

   	for(i in 1:N)
   	{
        r  <- r + rn;
		rn <- rn* (a + i) * (b + i) * z / (c + i) / (1 + i);
	}        
	
	return(r);
}

d_cor_dist_log<-function( x, rho.true, N.sample)
{
	N <- N.sample;
	rho <- rho.true;
	
	if(length(which(rho==1))>0) rho[which(rho==1)] <- 1 - 10^(-10);
	if(length(which(x==1))>0) x[which(x==1)] <- 1 - 10^(-10);
	
	x1 <- log(N-2) + ((N-1)/2)*log(1-rho^2) + ((N-4)/2)*log(1-x^2);
	x2 <- log( sqrt(2*pi)) + (N-1.5)* log(1-rho*x)

    gama.log <- lgamma(N-1) - lgamma(N-0.5);
    F  <- d_cor_dist_F( 0.5, 0.5, N - 0.5, (rho * x + 1) / 2, N );

	r.log <- gama.log + x1- x2 + log(F);
	
   	return( r.log )
}

clusterCorr_fast <- function( observed_cor_matrix, cluster_vector, return.grp=F ) 
{
	num_vertices = nrow(observed_cor_matrix)
	cluster_cor_mat <- array(NA, dim=c(num_vertices, num_vertices));

	max_clsuer <- max(cluster_vector);
	mat_sum  <- array(0, dim=c(max_clsuer, max_clsuer));
	mat_mnum <- array(0, dim=c(max_clsuer, max_clsuer));

	row.mat <- row(observed_cor_matrix)
	col.mat <- col(observed_cor_matrix)
	
	s1 <- rep(cluster_vector, num_vertices);
	s2 <- rep(cluster_vector, each=num_vertices);

	df <- data.frame( id=1:(num_vertices*num_vertices), r=s1, c=s2, val=c(observed_cor_matrix) );

	diag.idx <- which(c(diag(num_vertices)) ==1);
	df.diag <- df[diag.idx,]
	df <- df[ -diag.idx,,drop=F]
	
	#df2 <- sqldf("select r, c, avg(val) as avg from df group by r,c");
	if (!return.grp)
		df2 <- aggregate(df$val, by=list(r=df$r, c=df$c), FUN=mean)
	else
		df2 <- aggregate(df$val, by=list(r=df$r, c=df$c), FUN=mean)
	
	colnames(df2)<-c("r", "c", "avg")
	
	#df2[ which(df2$r==df2$c), 3] <- 0.99;
	if (!return.grp)
		df2[ which(df2$r!=df2$c), 3] <- df2[ which(df2$r!=df2$c), 3]/log(max_clsuer);
	
	#df0<- sqldf("select df.id, df.r, df.c, df2.avg from df left join df2 on df.r=df2.r and df.c=df2.c");
	df0 <- merge(x = df, y = df2, by = c("r", "c"), all.x=TRUE);
	df0 <- df0[,c("id", "r", "c", "avg"), drop=F]

	colnames(df.diag) <- c("id", "r", "c", "avg");
	df0 <- rbind( df0, df.diag);

	df0 <- df0[ order(df0[,1]), ]

	if(!return.grp)
	{
		cluster_cor_mat <- matrix(df0$avg, nrow=num_vertices);
		return(cluster_cor_mat)
	}
	else
	{
		
		grp.mat <- array( 1.0, dim=c(max(cluster_vector), max(cluster_vector)));		
		for(i in 1:NROW(df2))
			grp.mat[ df2[i,1],df2[i,2] ] <- df2[i,3];
		
		return(grp.mat);
	}
}


extend_hcluster<-function(hc, clusters, mat)
{
	hc.height0 <- hc$height; 
	hc.height  <- c();; 

	hc.merge0  <- hc$merge; 
	hc.merge   <- c();

	i.last.idx <- c();	
	for (i in 1:max(clusters))
	{
		i.clusters <- which(clusters==i);
		i.last<- -1*i.clusters[1];
		if(length(i.clusters )>1)
			for( j in 2:length(i.clusters ))
			{
				hc.merge <- rbind( hc.merge, c( -i.clusters[j], i.last) );
				hc.merge.idx <- NROW(hc.merge);
				i.last <- hc.merge.idx;
			
				hc.height <- c(hc.height, 0); 
			}
		
		i.last.idx <- c(i.last.idx, i.last);
	}
	
	hc.offset <- NROW(hc.merge);
	
	for(i in 1:NROW(hc.merge0))
	{
		hc.left <- hc.merge0[i,1];
		if( hc.left<0 )
			hc.left <- i.last.idx[-hc.left]
		else			
			hc.left <- hc.left + hc.offset;

		hc.right <- hc.merge0[i,2];
		if( hc.right<0 )
			hc.right <- i.last.idx[-hc.right]
		else			
			hc.right <- hc.right + hc.offset;

		hc.merge <- rbind( hc.merge, c( hc.left, hc.right) );
		hc.height <- c(hc.height, hc.height0[i]); 
	}

	hc.order0 <- hc$order;
	hc.order <- c();
	
	for(i in 1:length(hc.order0))
	{
		hc.order <- c(hc.order, which(clusters==hc.order0[i]) );
	}
	
	hc$height <- hc.height;
	hc$order <- hc.order;
	hc$merge <- hc.merge;
	
	return(hc);
}

tfbs_drawheatmapForClusters<-function( tfbs, mat, clusters ) 
{
	yb.sig.pal <- function(n, scale=10) 
	{
		ints<- c(0:(n-1))/(n-1) ## Linear scale from 0:1 x N values.
		ints<- 1/(1+exp(scale*(0.5-ints)))## Transfer to sigmoidal scale.

		b<- min(ints)
		m<- 2*b/(n-1)
		ints<- ints+(m*c(0:(n-1)) -b)## Transform by linear function to fill colors out to maxes.

		## Transfer to colorspace.
		# Yellow: 255, 255, 0
		# White: 255, 255, 255
		# Blue: 0, 0, 255

		YW <- ints[ints < 0.5] *2
		WB <- (ints[ints >= 0.5]-0.5) *2
		YW[YW<0] <- 0; WB[WB>1] <- 1
		c(rgb(1, 1, YW), rgb(1-WB, 1-WB, 1))
	}


	mat.grp <- clusterCorr_fast(mat, clusters, T );
	hc0 <- agnes(as.dist((1-mat.grp)^5), diss=TRUE)

	hc <- extend_hcluster( hc0, clusters, mat)
	hc1 <- as.dendrogram(hc)
	ord.hc1 <- order.dendrogram(hc1)

    wts <-  diag(mat.grp)[clusters[hc$order]]
	hc2 <- reorder(hc1, wts)
	ord.hc2 <- order.dendrogram(hc2)

	cuth <- clusters[ ord.hc2 ];
	pal100 <- c("#7E291B","#66E52C","#8F66F0","#58DBE8","#396526","#EAABC1","#E1C33C","#3E3668","#EB3F90","#C3E6A8","#E74618","#66A2E9","#3E7774","#DF9056","#3C2C21","#DF40D7","#6CEF92","#8C5A6B","#BC8AE2","#A03B99","#56AC2D","#389C6C","#E26B7E","#706B4B","#D2E374","#A0A560","#7B1C3E","#49F7DB","#C8C6E9","#414FA2","#95A590","#8A669A","#98A62F","#9E792C","#D69489","#547FE5","#DF6340","#849BAB","#E63A61","#8A386D","#DAE338","#263715","#BBDEE3","#3F1324","#A1E03B","#383544","#76C2E8","#794D38","#DD74E2","#D7BF8E","#E366B9","#894D19","#D57221","#D9B764","#B0303D","#D6BFBC","#757A28","#D991CC","#356344","#E3A22D","#223C36","#83B664","#D8E4C5","#AE9ED9","#37576D","#CF7398","#6BEAB2","#6C9266","#B33662","#4B340E","#57E65B","#E23635","#9B464B","#757089","#578BBA","#A6311B","#B2714E","#457F20","#4EA3B1","#B93286","#73D463","#531914","#6B78C2","#E07467","#8D786E","#515018","#361E40","#AA4FCC","#90D2C4","#71469B","#419C4C","#37558A","#B393B0","#9BE090","#856BDA","#66B593","#47CA84","#5D205B","#672F3F","#59D8C2")
	pal500 <- rep(pal100, 5)

	fill.col<- rep("white", length(ord.hc2));
	fill.col[ord.hc2] <- pal500[cuth];
	
	print( levelplot((mat)[ord.hc2, ord.hc2], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="",
		colorkey = list(space="left", labels=list(cex=1.5)), 
		legend = list(
          top = list(fun = dendrogramGrob,
          args = list(x = hc2, side = "top", #ord = ord.hc2, #lwd=2,
          size = 7, size.add = 0.5, 
          add = list(rect = list(col = "transparent", fill= fill.col )),
          type = "rectangle")))) );
}
