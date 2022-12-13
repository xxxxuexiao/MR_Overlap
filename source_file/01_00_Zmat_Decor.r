# remove impact of correlation due to sample overlap on Zscore matrix
# @zmat: a n*d matrix of original z-scores with d equals to the number of SNPs common to n gwas studies; 
# @cor.mat: a n*n correlation matrix 

library(powerplus)

decor <- function(zmat, cor.mat)
{
	zdecor = Matpow(cor.mat,-0.5) %*% zmat		
	dimnames(zdecor) <- dimnames(zmat)
	
	return(zdecor)	
}
