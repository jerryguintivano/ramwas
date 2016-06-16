rowSumsSq = function(mat) {
	output = double(NROW(mat));
	.Call("CrowSumsSq", mat, output);
	return(output);
}

colSumsSq = function(mat) {
	output = double(NROW(mat));
	.Call("CcolSumsSq", mat, output);
	return(output);
}
