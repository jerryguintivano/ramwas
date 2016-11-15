# Find "CG"s in the reference genome
getCpGsetCG = function( genome ){
    # library(Biostrings)
    cpgset = vector("list", length(genome))
    names(cpgset) = names(genome);
    for( i in seq_along(genome) ){ # i = length(genome)
        cpgset[[i]] = start(matchPattern("CG", genome[[i]], fixed=TRUE));
    }
    return(cpgset);
}

# Find any pair of letters which can become "CG"
# in reference genome with SNPs injected
getCpGsetALL = function( genome ){
    Cset = c("C","Y","S","M","B","H","V"); #,"N"
    Gset = c("G","R","S","K","B","D","V"); #,"N"

    Craw = sapply(Cset,charToRaw);
    Graw = sapply(Gset,charToRaw);

    Cint = logical(256); Cint[as.integer(Craw)] = TRUE;
    Gint = logical(256); Gint[as.integer(Graw)] = TRUE;

    starts_fun = function(chrom){
        which( Cint[as.integer(chrom[-length(chrom)])] &
               Gint[as.integer(chrom[-1])] )
    }

    cpgset = vector("list", length(genome))
    names(cpgset) = names(genome);
    for( i in seq_along(genome) ){ # i = length(genome)
        message("Processing ",names(genome)[i]);
        cpgset[[i]] = starts_fun(charToRaw(as.character(genome[[i]])));
    }
    return(cpgset);
}

# Create compressed FASTQ files
# in in-silico alignment experiment
insilicoFASTQ = function(con, gensequence, fraglength){
    # con=""; gensequence = "ABCDEFG"; fraglength=4;
    # con="D:/fastq.gz"; gensequence = "ABCDEFG"; fraglength=4;

    if (is.character(con)){
        if(nchar(con) > 0){
            if(grepl("\\.gz$",con)){
                con = gzfile(con, open = "wb")
            } else {
                con = file(con, open = "wb")
            }
            on.exit(close(con))
        } else {
            con = NULL;
        }
    }

    qual = charToRaw(paste0("\n+\n",paste(rep("A",fraglength),
                                          collapse = ""),
                            "\n"));

    sequence = as.character(gensequence);
    Encoding(sequence) = "bytes";
    y = charToRaw(sequence);
    rm(sequence);
    len = length(y);

    mat = NULL;

    step1 = 102400;
    mm = len - fraglength+1;
    nsteps = ceiling(mm/step1);
    for( part in seq_len(nsteps) ){ # part=1
        if(!is.null(con))
            message("step ", part, " of ", nsteps);
        fr = (part-1)*step1 + 1;
        to = min(part*step1, mm);
        if( NCOL(mat) != (to-fr+1) ){
            mat = vector("list",3*(to-fr+1));
            dim(mat) = c(3,to-fr+1);
            mat[3,] = list(qual);
        }
        mat[1,] = lapply(paste0("@",formatC(fr:to, width = 9, flag = "0"),"\n"),
                         charToRaw);
        mat[2,] = lapply(fr:to,
                         function(a){ y[a:(a+fraglength-1)] } );

        keep = (y[fr:to]!=0x4e) & (y[(fr:to)+fraglength-1]!=0x4e);
        if(any(keep)){
            if(is.null(con)){
                cat(rawToChar(unlist(mat[,keep])));
            } else {
                writeBin(con = con, object = unlist(mat[,keep]));
            }
        }
    }
    rm(part, step1, mm, nsteps, fr, to);
    return(invisible(TRUE));
}

# inject SNPs into a genome sequence
# with MAF filtering
injectSNPsMAF = function(gensequence, frqcount, MAF = 0.01){ #
    # http://droog.gs.washington.edu/parc/images/iupac.html
    {
        strACGT = array("", c(2,2,2,2));

        strACGT[1,1,1,1] = "*";

        strACGT[2,1,1,1] = "A";
        strACGT[1,2,1,1] = "C";
        strACGT[1,1,2,1] = "G";
        strACGT[1,1,1,2] = "T";

        strACGT[2,2,1,1] = "M";
        strACGT[2,1,2,1] = "R";
        strACGT[2,1,1,2] = "W";
        strACGT[1,2,2,1] = "S";
        strACGT[1,2,1,2] = "Y";
        strACGT[1,1,2,2] = "K";

        strACGT[2,2,2,1] = "V";
        strACGT[2,2,1,2] = "H";
        strACGT[2,1,2,2] = "D";
        strACGT[1,2,2,2] = "B";

        strACGT[2,2,2,2] = "N";

        rawACGT = sapply(strACGT, charToRaw)
        dim(rawACGT) = dim(strACGT)

        ACGT = 1:4
        ACGTnames = c("A","C","G","T")
        names(ACGT) = ACGTnames
    }
    # chrom = charToRaw( readRDS(filechr) );
    if( length(frqcount)==1 ){
        frqcount = readLines(con = frqcount);
    }
    if(grepl("^CHROM\t",frqcount[1]))
        frqcount = frqcount[-1];

    spt = strsplit(x = frqcount, split = "\t", fixed = TRUE);
    rm(frqcount);
    gc();

    pos = as.integer( sapply( spt, `[`, 2) );

    gensequence = charToRaw(as.character(gensequence));

    # options(warn=2)
    for( i in seq_along(spt)){ # i=17
        if( (i %% 10000) == 0 )
            message("Step ", i, " of", length(spt));
        tl = spt[[i]][-(1:4)]
        tlspt = strsplit(tl,":",TRUE);
        allele = sapply(tlspt,`[`,1);
        count = as.integer(sapply(tlspt,tail,1));
        names(count) = allele;
        countACGT = count[ACGTnames];
        countACGT[is.na(countACGT)] = 0;
        if(sum(countACGT) > 0 ){
            countACGT = countACGT / sum(countACGT);
            index = (countACGT >= MAF)+1L;
            dim(index) = c(1,4);
            gensequence[pos[i]] = rawACGT[index];
        }
    }
    return(rawToChar(gensequence));
}
