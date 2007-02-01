.packageName <- "ARES";

"read.genepop" <- function( filename = NULL )
{

	genepop <- vector("list", 1);

	# read the file into memory
	con <- file( description = filename, open = "rt" );
	lines <- readLines( con );
	close( con );
	fileLength = length(lines);

	#pop locations ignoring first 'title' line.
	popLocations <- grep("POP", lines[2:fileLength], ignore.case=TRUE);
	popLocations <- popLocations + 1; 

	#resize output
	genepop <- vector("list", length(popLocations));

  	#get title from first line
  	title <- lines[1];

	#get Loci Column names
  	lociNames <- lines[2:(popLocations[1]-1)];
	lociNames <- gsub("\t","",lociNames);
	
	#parsing of data between pops into genepop
	for( i in 1:length(popLocations) )
	{
		beginLine <- popLocations[i] + 1;
		endLine <- 0;

		group <- vector("list", length(lociNames));

		if( i == length(popLocations) )
		{	
			endLine <- length(lines) - 1; 
		}
		else
		{	
			endLine <- popLocations[i+1] - 1; 
		}
				
		#parse individual line
		for( line in lines[beginLine:endLine] )
		{
			
			#split id & alleles apart
			location <- unlist(strsplit( line, ","));
			id <- location[1];
			alleleGroup <- location[2];

			#split alleles apart on whitespace
			alleles <- unlist(strsplit( alleleGroup, '[[:space:]]+', extended=TRUE ))
			alleles <- alleles[2:length(alleles)];

			#cat(" ",length(alleles)," => ", alleles ,"\n");	

			if( length(lociNames) == length(alleles) )
			{

				for( j in 1:length(alleles) )
				{
					
					#only recognize data that isn't missing
					if( ( alleles[j] != "0000" ) && ( alleles[j] != "000000" ) )
					{
						digits <- unlist(strsplit( alleles[j], "") );
						
						if( length(digits) == 4 )
						{
							digits1 <- unlist(substring( alleles[j], 1, 2 ));
							digits2 <- unlist(substring( alleles[j], 3, 4 ));
							if( digits1 != "00" )
							{
								group[[j]] <- c( group[[j]], digits1 );
							}
							if( digits2 != "00" )
							{
								group[[j]] <- c( group[[j]], digits2 );
							}
						}
						else if( length(digits) == 6 )
						{
							digits1 <- unlist(substring( alleles[j], 1, 3 ));
							digits2 <- unlist(substring( alleles[j], 4, 6 ));
							if( digits1 != "000" )
							{
								group[[j]] <- c( group[[j]], digits1 );
							}
							if( digits2 != "000" )
							{
								group[[j]] <- c( group[[j]], digits2 );
							}
						}
						
						group[[j]] <- unique( group[[j]] );
						group[[j]] <- sort( group[[j]] );
						
						allAttributes <- list( name=lociNames[j], id=id );
						attributes( group[[j]] ) <- allAttributes;
					
					
					}
					
				}

			}

		}
		
		attributes(group) <- list( id=i, beginLine=beginLine, endLine=endLine );
		genepop[[i]] <- group;

	}

	contents <- list( title=title, filename=filename, lines=lines );
	attributes( genepop ) <- contents;	

	title <- genepop@title;
	lines <- genepop@lines;
	
	#parsing of data between pops into genepop
	for( i in 1:length(genepop) )
	{
		
		group <- genepop[[i]];
		id <- attr(group,"id");
		beginLine <- attr(group,"beginLine");
		endLine <- attr(group,"endLine");
		
		colCount <- endLine - beginLine + 1;
		col <- 0;
		rowCount <- 0;
		rowNames = list;

		for( j in 1:length(genepop[[i]]) )
		{
			locus <- genepop[[i]][[j]];
			rowCount <- rowCount + length(locus);
			locusName <- attr(locus,"name");
			for( k in 1:length(locus) )
			{
				name <- paste( list(locusName,"-",locus[k]), collapse="") ;
				rowNames <- c( rowNames, name);

			}
		}
		
		rowNames <- rowNames[2:length(rowNames)];
		output <- matrix(0,nrow=rowCount, ncol=colCount, byrow=TRUE );
		
		#parse individual line
		for( line in lines[beginLine:endLine] )
		{
			
			#split id & alleles apart
			location <- unlist(strsplit( line, ","));
			id <- location[1];
			alleleGroup <- location[2];
			col <- col + 1;

			#split alleles apart on whitespace
			alleles <- unlist(strsplit( alleleGroup, '[[:space:]]+', extended=TRUE ))
			alleles <- alleles[2:length(alleles)];

			#cat(" col( ",col," )\n");
			for( j in 1:length(alleles) )
			{
				#only recognize data that isn't missing
				if( ( alleles[j] != "0000" ) && ( alleles[j] != "000000" ) )
				{
					digits <- unlist(strsplit( alleles[j], "") );
					
					if( length(digits) == 4 )
					{
						digits1 <- unlist(substring( alleles[j], 1, 2 ));
						digits2 <- unlist(substring( alleles[j], 3, 4 ));
					}
					else if( length(digits) == 6 )
					{
						digits1 <- unlist(substring( alleles[j], 1, 3 ));
						digits2 <- unlist(substring( alleles[j], 4, 6 ));
					}
					
					existDigits1 <- grep(digits1, group[[j]]);
					existDigits2 <- grep(digits2, group[[j]]);
					
					if( length(existDigits1) > 0 )
					{
						matchName <- paste( list(attr(group[[j]],"name"),"-",digits1), collapse="") ;
						rowNum <- grep(matchName, rowNames);
						if( length(rowNum) > 0 )
						{
							output[rowNum,col] <- output[rowNum,col] + 1;
						}
					}
					
					if( length(existDigits2) > 0 )
					{
						matchName <- paste( list(attr(group[[j]],"name"),"-",digits2), collapse="") ;
						rowNum <- grep(matchName, rowNames);
						if( length(rowNum) > 0 )
						{
							output[rowNum,col] <- output[rowNum,col] + 1;
						}
					
					}
				
				}
				
			}
			
			attr(genepop[[i]], "output") <- output;

		}

	}


	# return structure
	return( genepop );
}

