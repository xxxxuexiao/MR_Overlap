**01_bal_InSIDE_overlap_simulate_beta.R**: Individual-level simulations under balanced pleiotropy and InSIDE assumption satisfied.

**02_merge_individual_betahat.r**: merge all betahat for 100 simulations

**03_MR_analysis.r**: run MR using 10 different methods.

**04_eval_MR_analysis.r**: calculate MSE, type 1 error and power over 100 simulations

**source_file/02_03_Fit_tmvtnorm.r**: estimate overlap correlation coefficient between X and Y by fitting trucated multivariate normal distribution

**source_file/01_00_Zmat_Decor.r**: decorrelate z-score matrix with overlap correlation estimated by fitting trucated multivariate normal distribution

**source_file/EBay_PRS_linear.R**: calculate corrected beta_tweddie and fdrxbetatweddie

**source_file/EBay_PRS_linear.SE.r**: estimate SE for corrected fdrxbetatweddie using different bootstrapping methods
