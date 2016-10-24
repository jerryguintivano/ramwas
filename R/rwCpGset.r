getCpGset = function( genome ){
	cpgset = vector('list', length(genome))
	names(cpgset) = names(genome);
	for( i in seq_along(genome) ){ # i = length(genome)
		cpgset[[i]] = start(matchPattern('CG', genome[[i]], fixed=TRUE));
	}
	return(cpgset);
}

insilicoFASTQ = function(con, gensequence, fraglength){
	# con=""; gensequence = "ABCDEFG"; fraglength=4;
	# con="D:/fastq.gz"; gensequence = "ABCDEFG"; fraglength=4;
	
	if (is.character(con)) {
		if(nchar(con) > 0) {
			if(grepl('\\.gz$',con)) {
				con = gzfile(con, open = 'wb')
			} else {
				con = file(con, open = 'wb')
			}
			on.exit(close(con))
		} else {
			con = NULL;
		}
	}
	
	qual = charToRaw(paste0('\n+\n',paste(rep('A',fraglength),collapse = ''),'\n'));

	sequence = as.character(gensequence);
	Encoding(sequence) = "bytes";
	y = charToRaw(sequence);
	rm(sequence);
	len = length(y);
	
	mat = NULL;
	
	step1 = 102400;
	mm = len - fraglength+1;
	nsteps = ceiling(mm/step1);
	for( part in seq_len(nsteps) ) { # part=1
		if(!is.null(con))
			message('step ', part, ' of ', nsteps);
		fr = (part-1)*step1 + 1;
		to = min(part*step1, mm);
		if( NCOL(mat) != (to-fr+1) ) {
			mat = vector('list',3*(to-fr+1));
			dim(mat) = c(3,to-fr+1);
			mat[3,] = list(qual);
		}
		mat[1,] = lapply(paste0('@',formatC(fr:to, width = 9, flag = "0"),'\n'), charToRaw);
		mat[2,] = lapply(fr:to, function(a) { y[a:(a+fraglength-1)] } );
		
		keep = (y[fr:to]!=0x4e) & (y[(fr:to)+fraglength-1]!=0x4e);
		if(any(keep)) {
			if(is.null(con)) {
				cat(rawToChar(unlist(mat[,keep])));
			} else {
				writeBin(con = con, object = unlist(mat[,keep]));
			}
		}
	}
	rm(part, step1, mm, nsteps, fr, to);
	return(invisible(TRUE));
}


