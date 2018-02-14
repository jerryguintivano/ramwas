# Data matrix access class
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
        initialize = function(param = NULL, getPCs = TRUE, lockfile = NULL){
            if( !is.null(param) ){
                .self$open(param = param, getPCs = getPCs, lockfile = lockfile);
                return(invisible(.self));
            }
            .self$fmdata = new("filematrix");
            .self$samplenames = character(0);
            .self$nsamples = 0;
            .self$ncpgs = 0;
            .self$ndatarows = 0;
            .self$rowsubset = NULL;
            .self$cvrtqr = NULL;
            return(invisible(.self));
        },
        open = function(param, getPCs = TRUE, lockfile = NULL){
            # Checks of parameters and files
            
            # Covariates defined
            if(is.null(param$covariates))
                stop(
                    "Covariates are not defined.\n",
                    "See \"filecovariates\" or \"covariates\" parameter.");
            
            # All covariates present
            cvrtset = match(
                        x = param$modelcovariates,
                        table = names(param$covariates), 
                        nomatch = 0L);
            if( any(cvrtset == 0L) )
                stop( 
                    "The \"modelcovariates\" lists unknown covariates: \n",
                    paste0(param$modelcovariates[head(which(cvrtset==0))], 
                        collapse = ", "));
            
            # Extract covariates
            cvrt = param$covariates[ cvrtset ];
            rm(cvrtset);
            
            if( any(sapply(lapply(cvrt, is.na), any)) )
                stop("Missing values are not allowed in the covariates.");

            
            # Sample names in covariates
            .self$samplenames = as.character(param$covariates[[1]]);
            .self$nsamples = length(samplenames);
        
            # Open data matrix
            .self$fmdata = fm.open( 
                    filenamebase = paste0(param$dircoveragenorm, "/Coverage"), 
                    readonly = TRUE,
                    lockfile = lockfile);
            fmsamples = rownames(fmdata);
            .self$ncpgs = ncol(fmdata);
            .self$ndatarows = nrow(fmdata);
            # nsamplesall = nrow(fmdata);
        
            # Match samples in covariates with those in coverage matrix
            .self$rowsubset = match(samplenames, fmsamples, nomatch = 0L);
            if( any(rowsubset == 0L) )
                stop( "Unknown samples in covariate file: ",
                    paste(samplenames[head(which(rowsubset==0))],
                        collapse = ", "));
        
            # if no reordering is required, set rowsubset=NULL
            if( length(samplenames) == length(fmsamples) ){
                if( all(rowsubset == seq_along(rowsubset)) ){
                    .self$rowsubset = NULL;
                }
            }
            
            # Get PCs
            if( getPCs & (param$modelPCs > 0) ){
                filename = paste0(param$dirpca, "/eigenvectors.bmat");
                if( !file.exists(filename) )
                    stop(   "File not found: ", filename, "\n",
                            "Cannot include PCs in the analysis.\n",
                            "Run PCA analysis first with ramwas4PCA().");
                # e = readRDS(filename);
                eigenvectors = fm.open(
                        filenamebase = paste0(param$dirpca, "/eigenvectors"),
                        readonly = TRUE);
                PCs = eigenvectors[, seq_len(param$modelPCs)]; # , drop=FALSE
                eigenvectors$close();
                
                if(!is.null( rowsubset ))
                    PCs = PCs[rowsubset,];
                cvrt = cbind(cvrt, PCs);
                rm(e);
            }
        
            .self$cvrtqr = orthonormalizeCovariates(
                                cvrt = cvrt,
                                modelhasconstant = param$modelhasconstant);
            
            return(invisible(.self));
        },
        close = function(){
            fmdata$close();
        },
        getDataRez = function(colset, resid = TRUE){
            # Get data
            x = fmdata[, colset];
            
            # Subset to active rows
            if( !is.null(rowsubset) )
                x = x[rowsubset, ];
            
            # NA's are coded with 255 in raw filematrices
            if( is.raw(x) ){
                # x = matrix(as.raw(0:255), 16, 16);
                storage.mode(x) = 'numeric';
                x[x == 255] = NA_real_;
            }
            
            # Impute missing values
            naset = is.na(x);
            if( any(naset) ){
                set = which(colSums(naset) > 0L);
                rm(naset);
                for( j in set ){ # j = set[1]
                    cl = x[,j];
                    mn = mean(cl, na.rm = TRUE);
                    if( is.na(mn) )
                        mn = 0;
                    x[is.na(cl), j] = mn;
                }
            } else {
                rm(naset);
            }

            # Orthogonalize w.r.t. covariates
            if( resid ){
                x = x - tcrossprod(cvrtqr, crossprod(x, cvrtqr));
            }
            
            return(x);
        }
    )
)
