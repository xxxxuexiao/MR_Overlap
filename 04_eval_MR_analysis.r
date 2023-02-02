eval_MR <- function(dat, true_theta){
	dat = dat[rowSums(is.na(dat))==0, , drop=FALSE]
	
	# get mse
	mse = mean((dat[,1]-true_theta)^2)
	
	# get bias
	bias = mean(dat[,1]) - true_theta
	
	# get variance
	variance = mean((dat[,1] - mean(dat[,1]))^2)
	
	# get type 1 error or power
	if(true_theta==0){
		FP = sum(dat[,2]<0.05)
		TN = sum(dat[,2]>0.05)
		alpha = FP/(FP+TN)
	} else{	
		TP = sum(dat[,2]<0.05)
		FN = sum(dat[,2]>0.05)
		alpha = TP/(TP+FN)
	}
	
	return(list(mse, alpha, bias, variance))
}

eval_MR_summary <- function(input_path, true_theta){
	load(input_path)

	# extract results using different mr methods
	IVW = est[, c(3+1, 3+33)]
	median = est[, c(3+2, 3+34)]
	mode = est[, c(3+3, 3+35)]
	PRESSO = est[, c(3+4, 3+36)]
	Robust = est[, c(3+5, 3+37)]
	Lasso = est[, c(3+6, 3+38)]
	egger = est[, c(3+7, 3+39)]
	conmix = est[, c(3+8, 3+40)]
	MRMix = est[, c(3+9, 3+41)]
	RAPS = est[, c(3+10, 3+42)]

	# cal. mse, true positive rate(power) or false positive rate(type 1 error) for each mr method
	res = lapply(list(IVW, median, mode, PRESSO, Robust, Lasso, egger, conmix, MRMix, RAPS), eval_MR, true_theta=true_theta)
	mse = do.call(cbind, lapply(res, `[[`, 1))
	alpha = do.call(cbind, lapply(res, `[[`, 2))
	bias = do.call(cbind, lapply(res, `[[`, 3))
	variance = do.call(cbind, lapply(res, `[[`, 4))
	
	header = c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS")
	colnames(mse) = paste0(header, "_mse")	
	colnames(alpha) = paste0(header, "_alpha")
	colnames(bias) = paste0(header, "_bias")
	colnames(variance) = paste0(header, "_variance")
	
	pattern = gsub(".rda", "", basename(input_path))
	output = cbind(pattern, mse, alpha, bias, variance)	
	return(output)
}

boxplot_dat_extract <- function(input_path)
{
	load(input_path)

	# extract results using different mr methods
	IVW = est[, c(3+1)]
	median = est[, c(3+2)]
	mode = est[, c(3+3)]
	PRESSO = est[, c(3+4)]
	Robust = est[, c(3+5)]
	Lasso = est[, c(3+6)]
	egger = est[, c(3+7)]
	conmix = est[, c(3+8)]
	MRMix = est[, c(3+9)]
	RAPS = est[, c(3+10)]
	
	dat = cbind(IVW, median, mode, PRESSO, Robust, Lasso, egger, conmix, MRMix, RAPS)
	colnames(dat) = c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS")
	
	return(dat)
}

# path for MR analysis datults
path_nodecor_mr = "/mnt/data/xue/Data/04_MR_Overlap/05_trial/01_bal_InSIDE/02_run_mr/01_nodecor_nowinner"
path_decor_mr = "/mnt/data/xue/Data/04_MR_Overlap/05_trial/01_bal_InSIDE/02_run_mr/02_decor_nowinner"
path_winner_mr = "/mnt/data/xue/Data/04_MR_Overlap/05_trial/01_bal_InSIDE/02_run_mr/03_nodecor_winner/02_tweedie"
path_winner_decor_mr = "/mnt/data/xue/Data/04_MR_Overlap/05_trial/01_bal_InSIDE/02_run_mr/04_decor_winner/02_tweedie"

# load mr results
all_thetas = c(0.2, 0)
all_prop_invalid = c(0.3)
outDir = "/mnt/data/xue/Data/04_MR_Overlap/05_trial/01_bal_InSIDE/02_run_mr/05_summary"

for (theta in all_thetas){
	for (prop in all_prop_invalid){
		path1 = paste0(path_nodecor_mr, "/X_Y25_theta", theta, "_thetaU0.5_N50000_prop_invalid", prop, ".rda")
		res1 = t(eval_MR_summary(path1, true_theta=theta))
		res1[1,1] = "no_decor_no_winner's_curse"
		dat1 = boxplot_dat_extract(path1)

		path2 = paste0(path_decor_mr, "/X_Y25_theta", theta, "_thetaU0.5_N50000_prop_invalid", prop, ".rda")
		res2 = t(eval_MR_summary(path2, true_theta=theta))
		res2[1,1] = "decor_no_winner's_curse"
		dat2 = boxplot_dat_extract(path2)

		# path3 = paste0(path_winner_mr, "/X_Y25_theta", theta, "_thetaU0.5_N50000_prop_invalid", prop, ".rda")
		# res3 = t(eval_MR_summary(path3, true_theta=theta))
		# res3[1,1] = "no_decor_winner's_curse"
		# dat3 = boxplot_dat_extract(path3)

		# path4 = paste0(path_winner_decor_mr, "/X_Y25_theta", theta, "_thetaU0.5_N50000_prop_invalid", prop, ".rda")
		# res4 = t(eval_MR_summary(path4, true_theta=theta))
		# res4[1,1] = "winner's_curse_decor"
		# dat4 = boxplot_dat_extract(path4)

		# res = cbind(res1, res2, res3, res4)
		res = cbind(res1, res2)
		write.table(res, file=paste0(outDir, "/X_Y25/X_Y25_theta", theta, "_thetaU0.5_N50000_prop_invalid", prop, ".summary.txt"), row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)
		
		header = c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS")
		# dat = lapply(header, function(method) { return(do.call(cbind, lapply(list(dat1, dat2, dat3, dat4), function(x) x[, method]))) })
		dat = lapply(header, function(method) { return(do.call(cbind, lapply(list(dat1, dat2), function(x) x[, method]))) })
		dat = do.call(cbind, dat)
		# colnames(dat) = paste0(rep(header, each=4), 1:4)
		colnames(dat) = paste0(rep(header, each=2), 1:2)
		write.table(dat, file=paste0(outDir, "/X_Y25/X_Y25_theta", theta, "_thetaU0.5_N50000_prop_invalid", prop, ".boxplot.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)		
	}
}


