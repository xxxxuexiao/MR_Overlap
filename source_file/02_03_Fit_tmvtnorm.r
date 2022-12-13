library(tmvtnorm)
library(data.table)
setDTthreads(20)
source("/mnt/data/xue/Data/02_COVID/01_data/00_Process_Sumstats.r")

# MLE fit Truncated Multivariate Normal Distribution(tmvtnorm) to get correlation under the null (i.e. correlation due to sample overlap).
fit.tmvtnorm <- function(gwasFilePaths, traits, merged_path=NA)
{ 
	common_snp = merge_sumstats(gwasFilePaths, traits, merged_path)
	zmat = as.matrix(common_snp[, c("SNP", grep("^Z_", colnames(common_snp), value=TRUE)), with=FALSE], rownames="SNP")
	zmat = zmat[rowSums(is.na(zmat))==0,]
	
 	# subset zmat with abs(Z)<1 to remove strongly associated snps, and remained snps are very weakly associated so that are treated as snps under the null (i.e. not associated with both phenotypes)
	indices = apply(zmat, 1, function(snp) all(abs(snp) < 1))
	zmat.fit = zmat[indices, ]
	
	# apply mle to fit truncated multivariate normal distribution
	lower = c(-1, -1)
	upper = c(1, 1)
	mle.fit = mle.tmvnorm(zmat.fit, lower=lower, upper=upper)
	
	# extract covariance from fitted model and cal. corrleation 
	# cor.xy = cov.xy / sqrt(cov.xx * cov.yy)
	cor_12 = as.numeric(coef(mle.fit)["sigma_1.2"] / sqrt(coef(mle.fit)["sigma_1.1"] * coef(mle.fit)["sigma_2.2"]))
	logL = as.numeric(logLik(mle.fit))
	
	return(c(cor_12, logL))
}
