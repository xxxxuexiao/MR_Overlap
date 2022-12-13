eval_MR <- function(dat, true_theta){
	dat = dat[rowSums(is.na(dat))==0, , drop=FALSE]
	
	mse = mean((dat[,1]-true_theta)^2)
	
	if(true_theta==0){
		FP = sum(dat[,2]<0.05)
		TN = sum(dat[,2]>0.05)
		alpha = FP/(FP+TN)
	} else{	
		TP = sum(dat[,2]<0.05)
		FN = sum(dat[,2]>0.05)
		alpha = TP/(TP+FN)
	}
	
	return(list(mse, alpha))
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
	colnames(mse) = paste0(c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS"), "_mse")
	colnames(alpha) = paste0(c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS"), "_alpha")
	
	pattern = gsub(".rda", "", basename(input_path))
	output = cbind(pattern, mse, alpha)	
	return(output)
}

# path for MR analysis datults
path_nodecor_mr = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/02_run_mr/01_nodecor_nowinner"
path_decor_mr = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/02_run_mr/02_decor_nowinner"
path_winner_mr = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/02_run_mr/03_nodecor_winner/02_tweedie"
path_winner_decor_mr = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/02_run_mr/04_decor_winner/02_tweedie"

# load mr results
all_thetas = c(0, 0.2, -0.2)
all_paths = c(path_nodecor_mr, path_decor_mr, path_winner_mr, path_winner_decor_mr)
all_prop_invalid = c(0.3, 0.5, 0.7)

for (theta in all_thetas){
	for (prop in all_prop_invalid){
		for (path in all_paths){			
			pattern = paste0("theta", theta, "_thetaU0.3_N50000_prop_invalid", prop)
			input_paths = list.files(path, pattern=pattern, full.names=TRUE)
			res = lapply(input_paths, eval_MR_summary, true_theta=theta)
			res = do.call(rbind.data.frame, res)
			write.table(res, file=paste0(path, "/", pattern, "_mr_mse_alpha.txt"), sep="\t", row.names=FALSE)
		}
	}	
}

outDir = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/02_run_mr/05_summary"

for (theta in all_thetas){
	for (prop in all_prop_invalid){
		path1 = paste0(path_nodecor_mr, "/X_Y75_theta", theta, "_thetaU0.3_N50000_prop_invalid", prop, ".rda")
		res1 = t(eval_MR_summary(path1, true_theta=theta))
		res1[1,1] = "no_decor_no_winner's_curse"

		path2 = paste0(path_decor_mr, "/X_Y75_theta", theta, "_thetaU0.3_N50000_prop_invalid", prop, ".rda")
		res2 = t(eval_MR_summary(path2, true_theta=theta))
		res2[1,1] = "decor_no_winner's_curse"

		path3 = paste0(path_winner_mr, "/X_Y75_theta", theta, "_thetaU0.3_N50000_prop_invalid", prop, ".rda")
		res3 = t(eval_MR_summary(path3, true_theta=theta))
		res3[1,1] = "no_decor_winner's_curse"

		path4 = paste0(path_winner_decor_mr, "/X_Y75_theta", theta, "_thetaU0.3_N50000_prop_invalid", prop, ".rda")
		res4 = t(eval_MR_summary(path4, true_theta=theta))
		res4[1,1] = "winner's_curse_decor"

		res = cbind(res1, res2, res3, res4)
		write.table(res, file=paste0(outDir, "/X_Y75/X_Y75_theta", theta, "_thetaU0.3_N50000_prop_invalid", prop, ".summary.txt"), row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)
	}
}