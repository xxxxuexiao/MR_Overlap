require(locfdr)
require(sfsmisc)

#***************************************************************************************************************************
# For computing empirical Bayes corrected estimates of regression coefficients for polygenic score testing/prediction
#***************************************************************************************************************************

#*******************************************************
# zval = list of z-statistics
# samp = sample size of the "base phenotype" ie the one you wish to construct PRS
# nulltype, df, bre: parameters for locfdr, please refer to the R package locfdr
#
# Output
# beta.corr.fdrXTweedie: corrected coefficient by tdr*Tweedie method (you can opt to include all variants in the polygenic score; ie choosing of p-value threhsolds can be avoided)
# beta.corr.Tweedie : corrected coefficient by Tweedie method
# beta.corr.fdr : corrected coefficient by tdr weighting method
#*******************************************************


# we assume standardized genotypes
EBay.PRS.linear <- function(zval, samp, nulltype=0, df=7, bre=120) {
  
  ztoVg.func <-function(Z,samp){
      Z^2 /(samp-2 + Z^2)
    }
  
  bias.corr.kernel<- function (zz,xout,bw="nrd0")
  {
    density.obj=density(zz,bw=bw) 
    fz.func = splinefun(density.obj$x, density.obj$y) 
    Psi.z = log(  fz.func(zz) /dnorm(zz)  )
    truez = D1ss(x= zz, y=Psi.z ,xout=xout)
    return(truez)
  }
  
  # Compute Tweedie formula adjusted z-values and beta
  zall.corr.ker =  bias.corr.kernel(zval,xout=zval)
  Vg.corr.ker = ztoVg.func(zall.corr.ker, samp  )
  beta.corr= sqrt( Vg.corr.ker )*sign(zval) 
  
  #Compute local fdr
  fdrobj = locfdr(zval, nulltype=nulltype, df=df, bre=bre, plot=0)
  fdr = fdrobj$fdr
  
  beta.standardized = ztoVg.func(zval, samp)
  beta.corr.fdr = (1-fdr)*beta.standardized
  beta.corr.fdrXTweedie = (1-fdr)*beta.corr
  return(list(beta.corr.fdrXTweedie = beta.corr.fdrXTweedie, 
              beta.corr.Tweedie = beta.corr,
              beta.corr.fdr = beta.corr.fdr))
}

##example
#zval =rnorm(1000)
#EBay.PRS.linear(zval, samp=20000)
