# Annotate findings with BioMaRt
ramwasAnnotateLocations = function(param, chr, pos){
    
    maxrequest = 500;
    
    # Sanity check
    if(any(pos > 1e9L))
        stop("Annotation error: chromosome positions must be <= 1e9")
    
    # BiomaRt likes 1-22,X,Y, not chr1-chr22,chrX,chrY
    if(is.character(chr)){
        nochr = gsub("^chr","",chr);
    } else {
        nochr = chr;
    }
    
    # Call biomaRt
    {
        # library(biomaRt)
        gene_ensembl = useMart(biomart=param$bimart,
                               host = param$bihost,
                               dataset=param$bidataset)
        
        if(length(nochr) > maxrequest){
            
            step1 = maxrequest;
            runto = length(nochr);
            nsteps = ceiling(runto/step1);
            collect = vector('list',nsteps)
            for( part in seq_len(nsteps) ) { # part = 1
                message('BiomaRt request: ', part, ' of ', nsteps, '\n');
                fr = (part-1)*step1 + 1;
                to = min(part*step1, runto);
                
                collect[[part]] = getBM(
                    mart = gene_ensembl,
                    attributes = unique(c(
                        param$biattributes,
                        "chromosome_name",
                        "start_position",
                        "end_position")),
                    filters = c(list(
                        chromosomal_region =paste0(
                            nochr[fr:to],":",
                            pos[fr:to] - param$biflank,":",
                            pos[fr:to] + param$biflank + 1L)),
                        param$bifilters));
            }
            rm(part, step1, runto, nsteps, fr, to);
            
            combine.data.frames = function(li) {
                rs = list();
                nms = names(li[[1]]);
                for( i in 1:length(nms) ) {
                    rs[[nms[i]]] = unlist(lapply(li,function(x){x[,nms[i]]}));
                }
                return(data.frame(rs));
            }
            
            bioresp = combine.data.frames(collect);
        } else {   
            bioresp = getBM(
                mart = gene_ensembl,
                attributes = unique(c(
                    param$biattributes,
                    "chromosome_name",
                    "start_position",
                    "end_position")),
                filters = c(list(
                    chromosomal_region =paste0(
                        nochr,":",
                        pos - param$biflank,":",
                        pos + param$biflank + 1L)),
                    param$bifilters));
        }
    }
    
    if(nrow(bioresp) == 0){
        rez = rep(list(rep("", length(chr))), length(param$biattributes))
        names(rez) = param$biattributes;
        return( data.frame(rez, stringsAsFactors = FALSE));
    }
    # Match Biomart response to chr/pos locations
    {
        # Transform chr/pos into single number coordinates
        chrunique = unique(nochr);
        
        tablepos = match(nochr, chrunique, nomatch = 0L) * 1e9 + pos;
        respchrn = match(bioresp$chromosome_name, chrunique, nomatch = 0L);
        resppos1 = respchrn * 1e9 + bioresp$start_position - param$biflank;
        resppos2 = respchrn * 1e9 + bioresp$end_position   + param$biflank;
        
        # Order CpGs by location
        ord = sort.list(tablepos);
        # tablepos[ord]
        
        ### Find CpGs which match each returned gene
        ### CpGs [fi1+1 : fi2] for each gene
        # offset by 1 for G in CpG, another 1 for strict inequality
        fi1 = findInterval( resppos1, tablepos[ord]+2)
        fi2 = findInterval( resppos2, tablepos[ord])
        
        ## Enumerate all CpG-gene pairs
        ## xx - CpG index (among sorted), factored for "split" call
        ## yy - Gene index (within biomart response)
        xx = unlist(lapply(which(fi1<fi2),
                           FUN=function(x){((fi1[x]+1):(fi2[x]))}));
        levels(xx) = paste0(seq_along(tablepos));
        class(xx) = "factor";
        yy = rep(seq_along(fi1), times = fi2-fi1)
        
        spl = split(yy, xx)
        
        result = vector("list",length(param$biattributes));
        names(result) = param$biattributes;
        
        for( attr in param$biattributes){ # attr = param$biattributes[1]
            z = bioresp[[attr]]
            result[[attr]] =
                sapply(spl, function(x){ paste0(z[x], collapse = "/") })
        }
        
        ord1 = seq_along(ord);
        ord1[ord] = ord1;
        
        # result2 = data.frame( chr, pos, lapply(result, `[`, ord1) )
        
        genes = data.frame(#chr = chr, 
            #pos = pos, 
            lapply(result, `[`, ord1), 
            stringsAsFactors = FALSE);	
    }
    return(genes);
}

ramwas6annotateTopFindings = function(param){
    # library(filematrix)
    param = parameterPreprocess(param);

    message("Working in: ", param$dirmwas);

    message("Loading MWAS results");
    mwas = fm.load( paste0(param$dirmwas, "/Stats_and_pvalues") );

    message("Loading CpG locations");
    cpgloc = fm.load(
        filenamebase = paste0(param$dircoveragenorm, "/CpG_locations") );
    chrnames = readLines(
        con = paste0(param$dircoveragenorm, "/CpG_chromosome_names.txt") );

    message("Finding top MWAS hits");
    keep = findBestNpvs(mwas[,3], param$toppvthreshold);
    # keep = which(mwas[,3] < param$toppvthreshold);
    ord = keep[sort.list(abs(mwas[keep,2]),decreasing = TRUE)];

    toptable = data.frame( chr = chrnames[cpgloc[ord,1]],
                           position =     cpgloc[ord,2],
                           tstat  = mwas[ord,2],
                           pvalue = mwas[ord,3],
                           qvalue = mwas[ord,4]);

    # saveRDS(file = paste0(param$dirmwas,"/Top_tests.rds"), object = toptable);

    if( !is.null(param$biattributes) && (nrow(toptable)>0L) ){
        message("Annotating top MWAS hits");
        bio = ramwasAnnotateLocations(param,
                                      chr = toptable$chr,
                                      pos = toptable$position);
        toptable = data.frame(toptable, bio);
    }

    message("Saving top MWAS hits");
    write.table(
        file = paste0(param$dirmwas,"/Top_tests.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        x = toptable
    );
    return(invisible(NULL));
}
