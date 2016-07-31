rowSumsSq = function(x) {
	stopifnot( is.numeric(x) );
	output = double(NROW(x));
	.Call("CrowSumsSq", x, output);
	return(output);
}

colSumsSq = function(x) {
	stopifnot( is.numeric(x) );
	output = double(NCOL(x));
	.Call("CcolSumsSq", x, output);
	return(output);
}
