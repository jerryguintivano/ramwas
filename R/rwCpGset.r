getCpGset = function( genome ){
	cpgset = vector('list', length(genome))
	names(cpgset) = names(genome);
	for( i in seq_along(genome) ){ # i = length(genome)
		cpgset[[i]] = start(matchPattern('CG', genome[[i]], fixed=TRUE));
	}
	return(cpgset);
}
