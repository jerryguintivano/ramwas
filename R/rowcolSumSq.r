# Same as rowSums(x^2), but faster with C/C++
rowSumsSq = function(x){
    stopifnot( is.numeric(x) );
    output = double(NROW(x));
    .Call("CrowSumsSq", x, output);
    return(output);
}

# Same as colSums(x^2), but faster with C/C++
colSumsSq = function(x){
    stopifnot( is.numeric(x) );
    output = double(NCOL(x));
    .Call("CcolSumsSq", x, output);
    return(output);
}
