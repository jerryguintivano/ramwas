if(FALSE){
    dr = "D:/temp"
    library(ramwas)
    param = ramwasParameters(
        dirproject = dr,
        dirbam = "bams",
        filebamlist = "bam_list.txt",
        filecpgset = "Simulated_chromosome.rds",
        cputhreads = 2,
        scoretag = "MAPQ",
        minscore = 4,
        minfragmentsize = 50,
        maxfragmentsize = 250,
        minavgcpgcoverage = 0.3,
        minnonzerosamples = 0.3,
        filecovariates = "covariates.txt",
        modelcovariates = NULL,
        modeloutcome = "age",
        modelPCs = 0,
        toppvthreshold = 1e-5,
        bihost = "grch37.ensembl.org",
        bimart = "ENSEMBL_MART_ENSEMBL",
        bidataset = "hsapiens_gene_ensembl",
        biattributes = c("hgnc_symbol","entrezgene","strand"),
        bifilters = list(with_hgnc_trans_name=TRUE),
        biflank = 0,
        cvnfolds = 5,
        mmalpha = 0,
        mmncpgs = c(5,10,50,100,500,1000,2000,3000)
    )
    
    ramwas4PCA(param);
    ramwas5MWAS(param);
    
    # param$modelcovariates = c('age','sex');
}
#
# param = parameterPreprocess(param);

setRefClass("rwDataClass",
	fields = list(
		fmdata = "filematrix",
		samplenames = "character",
		nsamples = "numeric",
		ncpgs = "numeric",
		ndatarows = "numeric",
		rowsubset = "ANY",
		cvrtqr = "ANY"
	),
	methods = list(
		initialize = function() {
			fmdata <<- new("filematrix");
		    samplenames <<- character(0);
		    nsamples <<- 0;
			ncpgs <<- 0;
			ndatarows <<- 0;
			rowsubset <<- NULL;
			cvrtqr <<- NULL;
			return(.self);
		},
		open = function( param, getPCs = TRUE, lockfile = NULL){
		    
		    # Checks of parameters and files
		    
		    # Covariates defined
		    if(is.null(param$covariates))
		        stop("Covariates are not defined.\n",
		             "See \"filecovariates\" or \"covariates\" parameter.");
		    
		    # All covariates present
		    cvrtset = match(param$modelcovariates, names(param$covariates), nomatch = 0L);
		    if( any(cvrtset == 0L) )
		        stop( "The \"modelcovariates\" lists unknown covariates: \n",
		              paste0(param$modelcovariates[head(which(cvrtset==0))], 
		                   collapse = ', '));
		    # Extract covariates
		    cvrt = param$covariates[ cvrtset ];
            rm(cvrtset);
		    
	        if( any(sapply(lapply(cvrt, is.na), any)) )
                stop("Missing values are not allowed in the covariates.");

		    
		    # Sample names in covariates
            samplenames <<- param$covariates[[1]];
            nsamples <<- length(samplenames);
        
            # Open data matrix
            fmdata <<- fm.open( 
                    filenamebase = paste0(param$dircoveragenorm, "/Coverage"), 
                    readonly = TRUE,
                    lockfile = lockfile);
            fmsamples = rownames(fmdata);
            ncpgs <<- ncol(fmdata);
            ndatarows <<- nrow(fmdata);
            # nsamplesall = nrow(fmdata);
        
            # Match samples in covariates with those in coverage matrix
            rowsubset <<- match(samplenames, fmsamples, nomatch = 0L);
            if( any(rowsubset == 0L) )
                stop( "Unknown samples in covariate file: ",
                    paste(samplenames[head(which(rowsubset==0))],
                        collapse = ', '));
        
            # if no reordering is required, set rowsubset=NULL
            if( length(samplenames) == length(fmsamples) ){
                if( all(rowsubset == seq_along(rowsubset)) ){
                    rowsubset <<- NULL;
                }
            }
            
            # Get PCs
            if( getPCs & (param$modelPCs > 0) ){
                filename = paste0(param$dirpca, "/eigen.rds");
                if( !file.exists(filename) )
                    stop(   "File not found: ", filename, "\n",
                            "Cannot include PCs in the analysis.\n",
                            "Run PCA analysis first with ramwas4PCA().");
                e = readRDS(filename);
                PCs = e$vectors[, seq_len(param$modelPCs), drop=FALSE];
                if(!is.null( rowsubset ))
                    PCs = PCs[rowsubset,];
                cvrt = cbind(cvrt, PCs);
                rm(e);
            }
        
            # Add a constant?
            if( param$modelhasconstant ){
                cvrt = c(const = list(rep(1, nrow(cvrt))), cvrt);
            } else {
                cvrt = cvrt;
            }
            
            # orthonormalize the covariates
            if( is.list(cvrt) ){
                isfactorset = sapply(cvrt, class) %in% c("character","factor");
                for( ind in seq_along(isfactorset) ){ # ind = 1
                    if( isfactorset[ind] ){
                        fctr = factor(cvrt[[ind]]);
                        if(nlevels(fctr) >= 2) {
                            cvrt[[ind]] = model.matrix(~fctr)[,-1];
                        } else {
                            cvrt[[ind]] = NULL;
                        }
                        rm(fctr);
                    } else {
                        # Kill pure zero covariates
                        if( all(cvrt[[ind]] == 0) )
                            cvrt[ind] = list(NULL);
                    }
                }
                cvrtmat = matrix(unlist(cvrt), nrow = nsamples);
            } else {
                cvrtmat = cvrt;
            }
            
            cvrtqr <<- qr.Q(qr(cvrtmat));
            
            
		},
		close = function(){
		    fmdata$close();
		},
		getDataRez = function(fr,to, resid = TRUE){
		    # Get data
		    x = fmdata[,fr:to];
		    
		    # Subset to active rows
		    if( !is.null(rowsubset) )
                x = x[rowsubset,];
		    
		    # Impute missing values
		    naset = is.na(x);
		    if( any(naset) ){
		        set = which(colSums(naset) > 0L);
		        for( j in set ) { # j = set[1]
		            cl = x[,j];
		            mn = mean(cl, na.rm = TRUE);
		            if( is.na(mn) )
		                mn = 0;
			        where1 = is.na( x[j, ] );
			        x[is.na(cl),j] = mn;
		        }
		    }
		    rm(naset);

		    # Orthogonalize w.r.t. covariates
		    if( resid ){
		        x = x - tcrossprod(cvrtqr, crossprod(x, cvrtqr));
		    }
		    
		    return(x);
		},
		getResults = function(){
		}
	)
)

# data = new("rwDataClass");
# data$open(param)
# data$cvrtqr
# data$getDataRez(1,10)

