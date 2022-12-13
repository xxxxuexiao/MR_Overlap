require(locfdr)
require(sfsmisc)
source("/mnt/data/xue/Data/04_MR_Overlap/00_script/EBay_PRS_linear.R")
 
EBay_PRS_linear.SE <- function(zall, method, totalN, repl=100) 
{
	# function to cal. Vg from a given z			                 
	ztoVg.func <-function(Z,samp=totalN){
		Z^2 /(samp-2 + Z^2)
	}

	# Modified Selection bias correction by Efron's empirical Bayes method       
	bias.corr.kernel<- function (zz,xout,bw="nrd0")
	{
		density.obj=density(zz,bw=bw)
		fz.func = splinefun(density.obj$x, density.obj$y)
		Psi.z = log(  fz.func(zz) /dnorm(zz)  )
		truez = D1ss(x= zz, y=Psi.z ,xout=xout)
		return(truez)
	}

	# parametric bootstrap with different resampling methods and weights
	# @method: resampling method, "paraboot", "fdrboot1" and "fdrboot2"
	boot.para <- function(listOfZ, repl, method){		   
		for(i in 1:repl) {
			z.sim = c() 
			for(j in 1:n){
				if(method == "paraboot"){
					z.sim[j] = rnorm(1, mean=listOfZ[j], sd=1)
				} else{
					dummy = runif(1, min=0,max=1); 
					if(dummy>fdr[j]){ z.sim[j]= rnorm(1,mean=listOfZ[j],sd=1) } else{ z.sim[j]=rnorm(1, mean=0,sd=1) }
				}				
			}									 				
			beta.cont.paraboot[[i]] = EBay.PRS.linear(zval=z.sim, samp=totalN)[[1]] # extract tdr*beta_twe
		}
		beta.cont.paraboot = do.call(rbind.data.frame, beta.cont.paraboot)
		boot.para.se = apply(beta.cont.paraboot, 2, sd)			
		return(boot.para.se)
	}

	Est.beta = EBay.PRS.linear(zval=zall, samp=totalN)[[1]] 
	
	beta.cont.paraboot = vector("list", length=repl)
	zall.corr.ker = bias.corr.kernel(zall,xout=zall)		
	fdrobj = locfdr(zall,bre=120,df=10,nulltype=2,plot=0)
	fdr = fdrobj$fdr
	n <- length(zall)
		
	# "parametric" bootstrap  
	# each sample New Z statistic = N(corrected z stastic,1)  
	if(method=="paraboot") {boot.se = boot.para(listOfZ = zall.corr.ker, repl=repl, method="paraboot")}
		
	# weighted "parametric" bootstrap 
	# each sample New Z statistic = N(observed Z stastic,1) with probability Pr(H1) 
	if(method=="fdrboot1") {boot.se = boot.para(listOfZ=zall, repl=repl, method="fdrtool1")}

	# weighted "parametric" bootstrap method 2 
	# each sample New Z statistic = N(corrected Z stastic,1) with probability Pr(H1)
	if(method=="fdrboot2") {boot.se = boot.para(listOfZ=zall.corr.ker, repl=repl, method="fdrboot2")}
	
	return(list(Est.beta=Est.beta, boot.se=boot.se)) 
}