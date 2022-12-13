# This program runs individual level simulation
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)

# path to save results
path_nodecor_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/01_nodecor_selected_IV"
path_decor_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/02_decor_selected_IV"
path_winner_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/03_winner_selected_IV/02_tweedie"
path_winner_decor_selectedIV = "/mnt/data/xue/Data/04_MR_Overlap/02_trial/01_bal_InSIDE/01_gwas_simulation/04_winner_decor_selected_IV/02_tweedie"


thetavec = c(0.2, 0, -0.2)
thetaU = 0.3
Nvec = c(5e4) # Nvec = c(5e4, 8e4, 1e5)
prop_invalid_vec = c(0.3, 0.5, 0.7)


merge_individual_betahat <- function(input, prefix)
{
	for (theta in thetavec){
		for (N in Nvec){
			for (prop_invalid in prop_invalid_vec){
				betahat_100sim = list()
				for (repgrp in 1:4){
					res = get(load(paste0(input,"/", prefix, "_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda")))
					betahat_100sim = c(betahat_100sim, res)
					rm(res)
				}
				print(length(betahat_100sim))
				save(betahat_100sim, file = paste0(input,"/100sim/betahat_100sim_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda"))
				rm(betahat_100sim)
			}
		}
	}
}

merge_individual_betahat(input=path_nodecor_selectedIV, prefix="betahat")
merge_individual_betahat(input=path_decor_selectedIV, prefix="betahat_decor")
merge_individual_betahat(input=path_winner_selectedIV, prefix="betahat_winner")
merge_individual_betahat(input=path_winner_decor_selectedIV, prefix="betahat_winner_decor")





