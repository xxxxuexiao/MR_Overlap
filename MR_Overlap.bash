# run different simulation paramaters for mr overlap
cd /mnt/data/xue/Data/04_MR_Overlap/05_trial/01_bal_InSIDE

# order in corresponding R vector
# thetavec=(1 2 3)
# thetaUvec=(1 2)
# Nvec=(1 2 3)
# prop_invalid_vec=(1 2 3)
# repgrp_vec=(1 2 3 4)

thetavec=(1 2 3)
thetaUvec=(1)
Nvec=(1)
prop_invalid_vec=(1 2 3)
repgrp_vec=(1 2 3 4)
idx_vec=(2 3 4 5)

# step1. gwas simulation
# 5 arguments represent "theta, thetaU, N, prop_invalid, repgrp"
for theta in "${thetavec[@]}"
do
  for thetaU in "${thetaUvec[@]}"
  do
    for N in "${Nvec[@]}"
    do
      for prop_invalid in "${prop_invalid_vec[@]}"
      do
        for repgrp in "${repgrp_vec[@]}"
        do
          echo "theta $theta thetaU $thetaU N $N prop_invalid $prop_invalid repgrp $repgrp start"
          Rscript 00_01_code/01_bal_InSIDE_overlap_simulate_beta.R $theta $thetaU $N $prop_invalid $repgrp > 00_02_log/theta"$theta"_thetaU"$thetaU"_N"$N"_prop_invalid"$prop_invalid"_repgrp"$repgrp".log 2>&1
		  echo "theta $theta thetaU $thetaU N $N prop_invalid $prop_invalid repgrp $repgrp done"
        done
      done
    done
  done
done

# step2. run mr
# 5 arguments represent "theta, thetaU, N, prop_invalid, idx"
for theta in "${thetavec[@]}"
do
  for thetaU in "${thetaUvec[@]}"
  do
    for N in "${Nvec[@]}"
    do
      for prop_invalid in "${prop_invalid_vec[@]}"
      do
        for idx in "${idx_vec[@]}"
        do
          echo "theta $theta thetaU $thetaU N $N prop_invalid $prop_invalid idx $idx start"
          Rscript 00_01_code/03_MR_analysis.r $theta $thetaU $N $prop_invalid $idx > 00_02_log/theta"$theta"_thetaU"$thetaU"_N"$N"_prop_invalid"$prop_invalid"_idx"$idx".log 2>&1
		  echo "theta $theta thetaU $thetaU N $N prop_invalid $prop_invalid idx $idx done"
        done
      done
    done
  done
done

