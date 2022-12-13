# Analyze summary statistics generated from individual level simulations
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)
library(MendelianRandomization)
library(MRMix)
library(mr.raps)
library(MRPRESSO)
library(penalized)

# path for betahat_100sim and corresponding mr results
path_nodecor_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/01_nodecor_selected_IV"
path_decor_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/02_decor_selected_IV"
path_winner_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/03_winner_selected_IV/02_tweedie"
path_winner_decor_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/04_winner_decor_selected_IV/02_tweedie"

path_nodecor_mr = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/02_run_mr/01_nodecor_nowinner"
path_decor_mr = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/02_run_mr/02_decor_nowinner"
path_winner_mr = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/02_run_mr/03_nodecor_winner/02_tweedie"
path_winner_decor_mr = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/02_run_mr/04_decor_winner/02_tweedie"


thetavec = c(0.2, 0, -0.2)
thetaUvec = c(0.3, 0.5)
Nvec = c(5e4, 8e4, 1e5)
prop_invalid_vec = c(0.3, 0.5, 0.7)

temp = as.integer(commandArgs(trailingOnly = TRUE))
# temp = c(1, 1, 1, 1, 2)
theta = thetavec[temp[1]] # True causal effect from X to Y
thetaU = thetaUx = thetaUvec[temp[2]] # Effect of the confounder on Y/X
N = Nvec[temp[3]] # Sample size for exposure X
prop_invalid = prop_invalid_vec[temp[4]] # Proportion of invalid IVs
idx = temp[5]# ith column for betahat_y of Y0, Y25, Y50 or Y75 for following MR analysis between X (idx=2,3,4,5)

pthr = 5e-8
NxNy_ratio = 1
M = 2e5

print(paste("N", N, "pthr", pthr, "theta", theta, "thetaU", thetaU, "prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio))

nx = N; ny = N/NxNy_ratio
est = matrix(NA, nrow = 100, ncol = 3+10*4+2)
mr_methods = c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS")
colnames(est) = c("numIV", "varX_expl","varY_expl", mr_methods, paste0(mr_methods,"_se"), paste0(mr_methods,"_time"), "PRESSO_pval", "conmix_CIcover0", paste0(mr_methods,"_P")) # varX_expl: variance of X explained by IVs; varY_expl: variance of Y explained by IVs

decor_vec = winner_vec = c(0, 1)
mynames = c("X", "Y0", "Y25", "Y50", "Y75")

for (decor in decor_vec){
	for (winner in winner_vec){
		## step1. load betahat_100sim
		if(decor==FALSE && winner==FALSE){
			load(paste0(path_nodecor_selectedIV,"/100sim/betahat_100sim_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda"))
			output_path = paste0(path_nodecor_mr, "/X_",mynames[idx],"_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda")	
			print(paste("run MR between X and", mynames[idx], "with no decor and no winner"))
		} else if(decor==TRUE && winner==FALSE){
			load(paste0(path_decor_selectedIV,"/100sim/betahat_100sim_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda"))
			output_path = paste0(path_decor_mr, "/X_",mynames[idx],"_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda")
			print(paste("run MR between X and", mynames[idx], "with decor but no winner"))
		} else if(decor==FALSE && winner==TRUE){
			load(paste0(path_winner_selectedIV,"/100sim/betahat_100sim_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda"))
			output_path = paste0(path_winner_mr, "/X_",mynames[idx],"_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda")
			print(paste("run MR between X and", mynames[idx], "with winner but no decor"))
		} else{
			load(paste0(path_winner_decor_selectedIV,"/100sim/betahat_100sim_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda"))
			output_path = paste0(path_winner_decor_mr, "/X_",mynames[idx],"_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda")
			print(paste("run MR between X and", mynames[idx], "with decor and winner"))
		}

		## step2. run mr
		for (repind in 1:100){
			set.seed(8967*repind)
			if(decor==FALSE) {numIV = nrow(betahat_100sim[[repind]])} else{numIV = nrow(betahat_100sim[[repind]][[idx-1]])}
			est[repind,1] = numIV
				
			if(numIV>2){
				if(decor==FALSE && winner==FALSE){		
					betahat_x.flt = as.vector(betahat_100sim[[repind]][,1])
					betahat_y.flt = as.vector(betahat_100sim[[repind]][,idx]) # idx=2,3,4,5 for Y0, Y25, Y50, Y75 
					bxse = rep(1/sqrt(nx), length(betahat_x.flt))
					byse = rep(1/sqrt(ny), length(betahat_y.flt))
				} else if(decor==TRUE && winner==FALSE){		
					betahat_x.flt = as.vector(betahat_100sim[[repind]][[idx-1]][,1]) # idx=1,2,3,4 for X_Y0, X_Y25, X_Y50, X_Y75
					betahat_y.flt = as.vector(betahat_100sim[[repind]][[idx-1]][,2])
					bxse = rep(1/sqrt(nx), length(betahat_x.flt))
					byse = rep(1/sqrt(ny), length(betahat_y.flt))
				} else if(decor==FALSE && winner==TRUE){		
					betahat_x.flt = as.vector(betahat_100sim[[repind]][,1])
					betahat_y.flt = as.vector(betahat_100sim[[repind]][,idx]) # idx=2,3,4,5 for Y0, Y25, Y50, Y75 
					bxse = as.vector(betahat_100sim[[repind]][,7])
					byse = as.vector(betahat_100sim[[repind]][,idx+6])
				} else{		
					betahat_x.flt = as.vector(betahat_100sim[[repind]][[idx-1]][,1]) # idx=1,2,3,4 for X_Y0, X_Y25, X_Y50, X_Y75
					betahat_y.flt = as.vector(betahat_100sim[[repind]][[idx-1]][,2])
					bxse = as.vector(betahat_100sim[[repind]][[idx-1]][,3])
					byse = as.vector(betahat_100sim[[repind]][[idx-1]][,4])
				}
				
				mr.obj = mr_input(bx = betahat_x.flt, bxse = bxse, by = betahat_y.flt, byse = byse)
				
				 # 1. IVW
				tryCatch({
					T0 = proc.time()[3]
					res = mr_ivw(mr.obj)
					T1 = proc.time()[3]
					est[repind,3+1] = res$Estimate
					est[repind,3+11] = res$StdError
					est[repind,3+21] = T1-T0
					est[repind,3+33] = res$Pvalue
					rm(res)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
				
				# 2. median
				tryCatch({
					T0 = proc.time()[3]
					res = mr_median(mr.obj)
					T1 = proc.time()[3]
					est[repind,3+2] = res$Estimate
					est[repind,3+12] = res$StdError
					est[repind,3+22] = T1-T0
					est[repind,3+34] = res$Pvalue
					rm(res)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
				
				# 3. mode
				tryCatch({
					T0 = proc.time()[3]
					res = mr_mbe(mr.obj)
					T1 = proc.time()[3]
					est[repind,3+3] = res$Estimate
					est[repind,3+13] = res$StdError
					est[repind,3+23] = T1-T0
					est[repind,3+35] = res$Pvalue
					rm(res)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
				
				# 4.MR-PRESSO
				tryCatch({
					presso.df = data.frame(bx = betahat_x.flt, by = betahat_y.flt, bxse = bxse, byse = byse)
					if (nx<=2e5){
						T0 = proc.time()[3]
						res = try(mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "byse", SdExposure = "bxse", 
											OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = presso.df, NbDistribution = 3000, SignifThreshold = 0.05))
						T1 = proc.time()[3]
						if (class(res)!="try-error"){
							if (!is.na(res$`Main MR results`[2,"Causal Estimate"]) & !is.na(res$`Main MR results`[2,"Sd"])){
								est[repind,3+4] = res$`Main MR results`[2,"Causal Estimate"]
								est[repind,3+14] = res$`Main MR results`[2,"Sd"]
								est[repind,34] = res$`Main MR results`[2,"P-value"]
								est[repind,3+36] = res$`Main MR results`[2,"P-value"]
							} else{
								est[repind,3+4] = res$`Main MR results`[1,"Causal Estimate"]
								est[repind,3+14] = res$`Main MR results`[1,"Sd"]
								est[repind,34] = res$`Main MR results`[1,"P-value"]
								est[repind,3+36] = res$`Main MR results`[1,"P-value"]
							}
							est[repind,3+24] = T1-T0
						} else{
							print(paste("Rep",repind,"numIV",numIV))
							print(res)
						}
						rm(res)
					}
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
				
				# 5. MR-Robust
				tryCatch({
					T0 = proc.time()[3]
					res = mr_ivw(mr.obj,"random", robust = TRUE)
					T1 = proc.time()[3]
					est[repind,3+5] = res$Estimate
					est[repind,3+15] = res$StdError
					est[repind,3+25] = T1-T0
					est[repind,3+37] = res$Pvalue
					rm(res)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
				
				# 6. MR-Lasso
				tryCatch({			
					T0 = proc.time()[3]
					res = mr_lasso(mr.obj)
					T1 = proc.time()[3]
					est[repind,3+6] = res$Estimate
					est[repind,3+16] = res$StdError
					est[repind,3+26] = T1-T0
					est[repind,3+38] = res$Pvalue
					rm(res)			
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
				
				# 7. Egger
				tryCatch({
					T0 = proc.time()[3]
					res = mr_egger(mr.obj)
					T1 = proc.time()[3]
					est[repind,3+7] = res$Estimate
					est[repind,3+17] = res$StdError.Est
					est[repind,3+27] = T1-T0
					est[repind,3+39] = res$Pvalue.Est
					rm(res)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
				
				# 8. contamination mixture
				tryCatch({
					T0 = proc.time()[3]
					res = mr_conmix(mr.obj)
					T1 = proc.time()[3]
					est[repind,3+8] = res$Estimate
					CIlength = res$CIUpper-res$CILower
					if (length(CIlength)>1) print(paste("Repind",repind,"conmix multimodal"))
					est[repind,3+18] = sum(CIlength)/1.96/2 ## Caution: this may be problematic
					est[repind,3+28] = T1-T0
					est[repind,35] = ifelse(sum((res$CILower<=0)&(res$CIUpper>=0))>0, 1, 0)
					est[repind,3+40] = res$Pvalue
					rm(res)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
				
				# 9. MRMix
				tryCatch({
					theta_temp_vec = seq(-0.5,0.5,by=0.01)
					T0 = proc.time()[3]
					res = MRMix(betahat_x.flt, betahat_y.flt, sx=1/sqrt(nx), sy=1/sqrt(ny), theta_temp_vec, pi_init = 0.6, sigma_init = 1e-5)
					T1 = proc.time()[3]
					est[repind,3+9] = res$theta
					est[repind,3+19] = res$SE_theta
					est[repind,3+29] = T1-T0
					est[repind,3+41] = res$pvalue_theta
					rm(res)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
				
				# 10. MR-RAPS
				tryCatch({
					T0 = proc.time()[3]
					res = mr.raps.overdispersed.robust(presso.df$bx, presso.df$by, presso.df$bxse, presso.df$byse,
													   loss.function = "huber", k = 1.345, initialization = c("l2"), 
													   suppress.warning = FALSE, diagnosis = FALSE, niter = 20, tol = .Machine$double.eps^0.5)
					T1 = proc.time()[3]
					est[repind,3+10] = res$beta.hat
					est[repind,3+20] = res$beta.se
					est[repind,3+30] = T1-T0
					est[repind,3+42] = res$beta.p.value
					rm(res)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})		
			} else{
				print(paste("Rep",repind,"numIV",numIV, "less than 2, skip"))
			}
				
			if (repind%%5==0){
				print(paste("Rep",repind,"numIV",numIV))
				save(est, file = output_path)
			}
		}
		save(est, file = output_path)	
	}
}




