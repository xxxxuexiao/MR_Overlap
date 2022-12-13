# simulations based on https://pubmed.ncbi.nlm.nih.gov/33393617/
# Individual-level simulations with balanced pleiotropy and InSIDE assumption satisfied. 

# step0: Generate indices of causal SNPs for four component
# step1: Simulate direct effect of G on X(gamma) or on Y(alpha) under M1
# step2: Simulate X, Y and overlapped Y across nbatch
		 # step2a. simulate X 
		 # step2b. simulate Y
		 # step2c. simulate overlapped Y
# step3: Filter the SNPs that reach genome-wide significance in the study associated with X. 
# step4: decorrelate with overlap correlation from mle fit of tmvtnorm
		 # step4a. convert betahat to z-score and decorrelate zmat
		 # step4b. filter IVs
# step5. only correct for winner's curse
		 # step5a. estimate corrected beta and corresponding se
			# (1) tweedie's formula corrected beta with unchanged SE(=sqrt(1/N))
			# (2) fdrxtweddie corrected beta with parametric bootstrap estimated SE
			# (3) fdrxtweddie corrected beta with fdrboot1 estimated SE
			# (4) fdrxtweddie corrected beta with fdrboot2 estimated SE
		 # step5b. filter IVs
# step6. first correct for winner's curse then decorrelate
	# (1)tweddie corrected beta + unchanged SE + decor
	# (2)tdr*tweddie corrected beta + paraboot SE + decor
	# (3)tdr*tweddie corrected beta + fdrboot1 SE + decor
	# (4)tdr*tweddie corrected beta + fdrboot2 SE + decor	
# Repeat simulation steps above 100 times


# simulation code start from here
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)
library(powerplus)
setDTthreads(20)
source("/mnt/data/xue/Data/02_COVID/04_script/02_03_Fit_tmvtnorm.r") # calculate overlap correlation between two traits
source("/mnt/data/xue/Data/02_COVID/03_cofdr/01_00_Zmat_Decor.r") # calculate decorrelated z-score matrix between two traits
source("/mnt/data/xue/Data/04_MR_Overlap/00_script/EBay_PRS_linear.R") # estimate corrected beta
source("/mnt/data/xue/Data/04_MR_Overlap/00_script/EBay_PRS_linear.SE.r") # estimate SE with bootstrap

# path to save results
path_raw_simualted_gwas = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/00_raw_simulated_gwas"
path_nodecor_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/01_nodecor_selected_IV"
path_decor_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/02_decor_selected_IV"
path_winner_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/03_winner_selected_IV"
path_winner_decor_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/04_winner_decor_selected_IV"


thetavec = c(0.2, 0, -0.2)
thetaUvec = c(0.3, 0.5)
Nvec = c(5e4, 8e4, 1e5)
prop_invalid_vec = c(0.3, 0.5, 0.7)

temp = as.integer(commandArgs(trailingOnly = TRUE)) # length=5
theta = thetavec[temp[1]] # True causal effect from X to Y
thetaU = thetaUx = thetaUvec[temp[2]] # Effect of the confounder on Y/X
N = Nvec[temp[3]] # Sample size for exposure X
prop_invalid = prop_invalid_vec[temp[4]] # Proportion of invalid IVs, pi2/(p1+i2)
repgrp = temp[5] # Simulation group number (1:8) - 200 simulations are partitioned into 8 groups for parallelization.

pthr = 5e-8 # p-value threshold for instrument selection
NxNy_ratio = 1 # Ratio of sample sizes for X and Y
M = 2e5 # Total number of independent SNPs representing the common variants between X and Y in the genome

# Model parameters for effect size distribution
pi1=0.02*(1-prop_invalid); pi3=0.01
pi2=0.02*prop_invalid; 
sigma2x = sigma2y = 5e-5; sigma2u = 0; sigma2x_td = sigma2y_td = (5e-5)-thetaU*thetaUx*sigma2u #if sigma2u!=0 then model M1 turns to M2

print(paste("N", N, "pthr", pthr, "pi1", pi1, "theta", theta, "thetaU", thetaU, "NxNy_ratio", NxNy_ratio))

# Sample size for X and Y
nx = N; ny = N/NxNy_ratio

# proportion of overlapping participants between X and Y
sample_overlap_prop = c(0.25, 0.5, 0.75)

# Due to large memory requirement of simulating genetic data for all the subjects at once, we simulate a batch of 500 subjects at a time, analyze each batch separately, and use meta-analysis to get final estimate.
batch_size = 500 
betahat_all = betahat_decor_all = vector("list", length = 25)
betahat_winner_all = betahat_winner_decor_all = vector("list", length = 25)
betahat_winner_all_fdrxtweddie_paraboot = betahat_winner_decor_all_fdrxtweddie_paraboot = vector("list", length = 25)
betahat_winner_all_fdrxtweddie_fdrboot1 = betahat_winner_decor_all_fdrxtweddie_fdrboot1 = vector("list", length = 25)
betahat_winner_all_fdrxtweddie_fdrboot2 = betahat_winner_decor_all_fdrxtweddie_fdrboot1 = vector("list", length = 25)

for (repind in 1:25){
	t1 = Sys.time()
    set.seed(89781*(repind+25*(repgrp-1)))
    
    # step0: Generate indices of causal SNPs for four component
    ind1 = sample(M, round(M*pi1))
    causalsnps = ind1
    ind2 = sample(setdiff(1:M,causalsnps), round(M*pi2))
    causalsnps = c(causalsnps,ind2)
    ind3 = sample(setdiff(1:M,causalsnps), round(M*pi3))
    causalsnps = c(causalsnps,ind3)

	
    # step1: Simulate direct effect of G on X(gamma) or on Y(alpha) under M1
    gamma = phi = alpha = rep(0,M)    
    gamma[ind1] = rnorm(length(ind1), mean = 0, sd = sqrt(sigma2x))
    gamma[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2x_td))
    alpha[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2y_td))
    alpha[ind3] = rnorm(length(ind3), mean = 0, sd = sqrt(sigma2y))
   

	# step2: Simulate X, Y and overlapped Y across nbatch
	nbatch = N/batch_size
	betahat_x = betahat_y = rep(0,M)
	betahat_y_overlap = vector("list", length=length(sample_overlap_prop)); betahat_y_overlap = lapply(betahat_y_overlap, function(x) x=rep(0,M))
	simu_gwas <- function(simu_X=TRUE, simu_Y=FALSE)
	{
		G = matrix(rbinom(batch_size*M, size=2, prob=0.3), ncol = M) # G of M snps with maf=0.3 for batch_size people
		G = (G-2*0.3)/sqrt(2*0.3*0.7) # standardized G with mean 0 and variance 1
		U = G %*% phi + rnorm(batch_size, mean = 0, sd = sqrt(1-M*pi2*sigma2u)) # simulate U with equation(1)
		X = G %*% gamma + thetaUx*U + rnorm(batch_size, mean=0, sd = sqrt(1-thetaUx^2-M*(pi1*sigma2x+pi2*sigma2x_td))) # simulate X with equation(2)   
		if(simu_X==TRUE){
			return(list(G=G, X=X))
		} else if(simu_Y==TRUE){
			Y = G %*% alpha + theta*X + thetaU*U + rnorm(batch_size, mean=0, sd = sqrt(1-theta^2-thetaU^2-2*theta*thetaU*thetaUx-M*(pi3*sigma2y+pi2*sigma2y_td))) # simulat Y with equation(3)
			return(list(G=G, Y=Y))
		} else{
			return(NULL)
		}
	}
	for (batch_ind in 1:nbatch){
		if (batch_ind%%5==0) print(batch_ind)
		
		# step2a. simulate X
		res_X = simu_gwas(simu_X=TRUE, simu_Y=FALSE)
		betahat_x = betahat_x + t(res_X$G)%*%res_X$X/batch_size # fit linear regression of G on X	

		# step2b. simulate Y
		res_Y = simu_gwas(simu_X=FALSE, simu_Y=TRUE)
		betahat_y = betahat_y + t(res_Y$G)%*%res_Y$Y/batch_size # fit linear regression of G on Y	
		
		# step2c. simulate overlapped Y
		overlap_ind = lapply(sample_overlap_prop, function(prop) sample(batch_size, round(batch_size*prop)))
		simu_Y_overlap <- function(ind, G_x=res_X$G, X=res_X$X, G_y=res_Y$G, Y=res_Y$Y)
		{
			G_y[overlap_ind[[ind]], ] <- G_x[overlap_ind[[ind]], ]
			Y[overlap_ind[[ind]], ] <- X[overlap_ind[[ind]], ]
						
			return(t(G_y)%*%Y/batch_size) # fit linear regression of overlapped_G on overlapped_Y
		}		
		for(i in 1:length(sample_overlap_prop)){			
			res_overlap = simu_Y_overlap(i)
			betahat_y_overlap[[i]] = betahat_y_overlap[[i]] + res_overlap
		}
		rm(res_X, res_Y)
	}
	betahat_x = betahat_x/nbatch
	betahat_y = cbind(betahat_y/nbatch, do.call(cbind, betahat_y_overlap)/nbatch, betahat_x)
	fwrite(betahat_y, file=paste0(path_raw_simualted_gwas,"/betahat_y_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repind",repind,"_repgrp",repgrp, ".txt.gz"))
 
    # step3: Filter the no-decorrelated SNPs that reach genome-wide significance in the study associated with X. 
    ind_filter = which(2*pnorm(-sqrt(nx)*abs(betahat_x))<pthr) 
	betahat.flt = cbind(betahat_x[ind_filter, , drop=FALSE], betahat_y[ind_filter, , drop=FALSE])	
	print(paste(c("rep", repind, "no decorrelation with", nrow(betahat.flt), "IVs including:", ind_filter), collapse=" "))
	betahat_all[[repind]] = betahat.flt
	save(betahat_all, file = paste0(path_nodecor_selectedIV,"/betahat_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))
	
	
	# step4. decorrelate with overlap correlation from mle fit of tmvtnorm
	common_snp = as.data.table(cbind(betahat_x*sqrt(nx), betahat_y*sqrt(ny)))	
	common_snp_se = matrix(sqrt(1/N), ncol=ncol(common_snp), nrow=nrow(common_snp))
	colnames(common_snp_se) = paste0("SE_", c("X", "Y0", "Y25", "Y50", "Y75", "Y100"))
	colnames(common_snp) = paste0("Z_", c("X", "Y0", "Y25", "Y50", "Y75", "Y100"))
	common_snp$SNP = paste0("snp", 1:nrow(common_snp))
	decorIV_filter <- function(t2, t1="X", merged_path_zval, merged_path_se){		
		## step4a. decorrelate
		traits = c(t1, t2)
		overlap_cor = fit.tmvtnorm(gwasFilePaths=NULL, traits=NULL, merged_path=merged_path_zval[, c("SNP", paste0("Z_",traits)), with=FALSE])[1] # calculate overlap correlation
		overlap_cor = matrix(c(1,overlap_cor,overlap_cor,1), nrow=2,ncol=2,byrow=TRUE)		
		zmat = as.matrix(merged_path_zval[, paste0("Z_",traits), with=FALSE])
		print(paste("decorrelate", t1, "and", t2,  "with overlap_cor"))
		print(overlap_cor)
		zmat_decor = t(decor(t(zmat), cor.mat=overlap_cor)) # calculate decorrelated z-score matrix 
		se_decor = merged_path_se[, paste0("SE_",traits)]
		betahat_decor = zmat_decor * se_decor # calculate decorrelated beta
		
		## step4b. filter IV based on decorrelated sumstats
		ind_filter_decor = which(2*pnorm(-abs(zmat_decor[,"Z_X"]))<pthr)
		betahat_decor.flt = betahat_decor[ind_filter_decor, , drop=FALSE] # drop=FALSE: treat a single row of matrix in R as a matrix object
		se_decor.flt = se_decor[ind_filter_decor, , drop=FALSE]
		colnames(betahat_decor.flt) = paste0("BETA_", traits)
		colnames(se_decor.flt) = paste0("SE_", traits)
		betahat_decor.flt = cbind(betahat_decor.flt, se_decor.flt)
		
		print(paste(c(nrow(betahat_decor.flt), "IVs including:", ind_filter_decor), collapse=" "))
		return(betahat_decor.flt)	
	}
	betahat_decor_all[[repind]] = lapply(c("Y0", "Y25", "Y50", "Y75", "Y100"), decorIV_filter, t1="X", merged_path_zval=common_snp, merged_path_se=common_snp_se)
	save(betahat_decor_all, file = paste0(path_decor_selectedIV, "/betahat_decor_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))


	# step5. correct for winner's curse
	# @method: beta correction method. "tweedie" or "fdrxtweddie"
	# @se.method: se estimation method for fdrxtweddie. "paraboot", "fdrboot1" or "fdrboot2"
	winnerIV_filter <- function(merged_path_zval, method="tweddie", se.method=NULL){
		## step5a. estimate corrected beta and corresponding se
		if(method == "tweddie"){
			res_winner = apply(merged_path_zval, 2, EBay.PRS.linear, samp=N)
			betahat_winner = do.call(cbind, lapply(res_winner, `[[`, 2))
			se_winner = matrix(sqrt(1/N), ncol=ncol(betahat_winner), nrow=nrow(betahat_winner)) # unchanged SE sqrt(1/N) for tweddie's formula corrected beta
		} else{
			res_winner = apply(merged_path_zval, 2, EBay_PRS_linear.SE, totalN=N, repl=100, method=se.method)
			betahat_winner = do.call(cbind, lapply(res_winner, `[[`, 1))
			se_winner = do.call(cbind, lapply(res_winner, `[[`, 2)) # bootstrap SE for fdrxtweddie corrected beta
		}
		z_winner = betahat_winner / se_winner
		colnames(betahat_winner) = c("BETA_X", "BETA_Y0", "BETA_Y25", "BETA_Y50", "BETA_Y75", "BETA_Y100")
		colnames(se_winner) = c("SE_X", "SE_Y0", "SE_Y25", "SE_Y50", "SE_Y75", "SE_Y100")
		
		## step5b. filter IVs
		ind_filter_winner = which(2*pnorm(-abs(z_winner[, "Z_X"]))<pthr)
		betahat_winner.flt = betahat_winner[ind_filter_winner, , drop=FALSE]
		se_winner.flt = se_winner[ind_filter_winner, , drop=FALSE]	
		betahat_winner.flt = cbind(betahat_winner.flt, se_winner.flt)
		
		print(paste(c("rep", repind, "using", method, "for winner correction with", nrow(betahat_winner.flt), "IVs including:", ind_filter_winner), collapse=" "))
		return(list(betahat_winner.flt, betahat_winner, se_winner))
	}
	## (1)tweddie corrected beta + unchanged SE
	winner_res = winnerIV_filter(merged_path_zval=common_snp[, 1:6], method="tweddie", se.method=NULL)
	betahat_winner_all[[repind]] = winner_res[[1]]
	save(betahat_winner_all, file = paste0(path_winner_selectedIV, "/02_tweedie/betahat_winner_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))	
	
	## (2)tdr*tweddie corrected beta + paraboot SE
	winner_res_fdrxtweddie_paraboot = winnerIV_filter(merged_path_zval=common_snp[, 1:6], method="fdrxtweddie", se.method="paraboot")
	betahat_winner_all_fdrxtweddie_paraboot[[repind]] = winner_res_fdrxtweddie_paraboot[[1]]
	save(betahat_winner_all_fdrxtweddie_paraboot, file = paste0(path_winner_selectedIV, "/01_fdrxtweedie/01_paraboot/betahat_winner_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))
	
	## (3)tdr*tweddie corrected beta + fdrboot1 SE
	winner_res_fdrxtweddie_fdrboot1 = winnerIV_filter(merged_path_zval=common_snp[, 1:6], method="fdrxtweddie", se.method="fdrboot1")
	betahat_winner_all_fdrxtweddie_fdrboot1[[repind]] = winner_res_fdrxtweddie_fdrboot1[[1]]
	save(betahat_winner_all_fdrxtweddie_fdrboot1, file = paste0(path_winner_selectedIV, "/01_fdrxtweedie/02_fdrboot1/betahat_winner_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))	
	
	## (4)tdr*tweddie corrected beta + fdrboot2 SE
	winner_res_fdrxtweddie_fdrboot2 = winnerIV_filter(merged_path_zval=common_snp[, 1:6], method="fdrxtweddie", se.method="fdrboot2")
	betahat_winner_all_fdrxtweddie_fdrboot2[[repind]] = winner_res_fdrxtweddie_fdrboot2[[1]]
	save(betahat_winner_all_fdrxtweddie_fdrboot2, file = paste0(path_winner_selectedIV, "/01_fdrxtweedie/03_fdrboot2/betahat_winner_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))	
	
		
	# step6. first correct for winner's curse then decorrelate
	# (1)tweddie corrected beta + unchanged SE + decor
	common_snp_winner = as.data.table(winner_res[[2]] / winner_res[[3]])
	common_snp_se_winner = winner_res[[3]]
	colnames(common_snp_winner) = paste0("Z_", c("X", "Y0", "Y25", "Y50", "Y75", "Y100"))
	common_snp_winner$SNP = paste0("snp", 1:nrow(common_snp_winner))		
	betahat_winner_decor_all[[repind]] = lapply(c("Y0", "Y25", "Y50", "Y75", "Y100"), decorIV_filter, t1="X", merged_path_zval=common_snp_winner, merged_path_se=common_snp_se_winner)
	save(betahat_winner_decor_all, file = paste0(path_winner_decor_selectedIV, "/02_tweedie/betahat_winner_decor_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))
	
	# (2)tdr*tweddie corrected beta + paraboot SE + decor
	common_snp_winner = as.data.table(winner_res_fdrxtweddie_paraboot[[2]] / winner_res_fdrxtweddie_paraboot[[3]])
	common_snp_se_winner = winner_res_fdrxtweddie_paraboot[[3]]
	colnames(common_snp_winner) = paste0("Z_", c("X", "Y0", "Y25", "Y50", "Y75", "Y100"))
	common_snp_winner$SNP = paste0("snp", 1:nrow(common_snp_winner))		
	betahat_winner_decor_all_fdrxtweddie_paraboot[[repind]] = lapply(c("Y0", "Y25", "Y50", "Y75", "Y100"), decorIV_filter, t1="X", merged_path_zval=common_snp_winner, merged_path_se=common_snp_se_winner)
	save(betahat_winner_decor_all_fdrxtweddie_paraboot, file = paste0(path_winner_decor_selectedIV, "/01_fdrxtweedie/01_paraboot/betahat_winner_decor_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))
	
	# (3)tdr*tweddie corrected beta + fdrboot1 SE + decor
	common_snp_winner = as.data.table(winner_res_fdrxtweddie_fdrboot1[[2]] / winner_res_fdrxtweddie_fdrboot1[[3]])
	common_snp_se_winner = winner_res_fdrxtweddie_fdrboot1[[3]]
	colnames(common_snp_winner) = paste0("Z_", c("X", "Y0", "Y25", "Y50", "Y75", "Y100"))
	common_snp_winner$SNP = paste0("snp", 1:nrow(common_snp_winner))		
	betahat_winner_decor_all_fdrxtweddie_fdrboot1[[repind]] = lapply(c("Y0", "Y25", "Y50", "Y75", "Y100"), decorIV_filter, t1="X", merged_path_zval=common_snp_winner, merged_path_se=common_snp_se_winner)
	save(betahat_winner_decor_all_fdrxtweddie_fdrboot1, file = paste0(path_winner_decor_selectedIV, "/01_fdrxtweedie/02_fdrboot1/betahat_winner_decor_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))
	
	# (4)tdr*tweddie corrected beta + fdrboot2 SE + decor
	common_snp_winner = as.data.table(winner_res_fdrxtweddie_fdrboot2[[2]] / winner_res_fdrxtweddie_fdrboot2[[3]])
	common_snp_se_winner = winner_res_fdrxtweddie_fdrboot2[[3]]
	colnames(common_snp_winner) = paste0("Z_", c("X", "Y0", "Y25", "Y50", "Y75", "Y100"))
	common_snp_winner$SNP = paste0("snp", 1:nrow(common_snp_winner))		
	betahat_winner_decor_all_fdrxtweddie_fdrboot2[[repind]] = lapply(c("Y0", "Y25", "Y50", "Y75", "Y100"), decorIV_filter, t1="X", merged_path_zval=common_snp_winner, merged_path_se=common_snp_se_winner)
	save(betahat_winner_decor_all_fdrxtweddie_fdrboot2, file = paste0(path_winner_decor_selectedIV, "/01_fdrxtweedie/03_fdrboot2/betahat_winner_decor_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))
	
	
	t2 = Sys.time()
	t.dur = difftime(t2, t1, units="hours")
	print(paste("rep", repind, "done using", t.dur, "hours"))
}