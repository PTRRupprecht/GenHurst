####################################
# Calculates the generalized Hurst exponent H(q) from the scaling 
# of the renormalized q-moments of the distribution 
#
#       <|x(t+r)-x(t)|^q>/<x(t)^q> ~ r^[qH(q)]
#
####################################
# H = genhurst(S,q)
# S is Tx1 data series (T>50 recommended)
# calculates H, specifies the exponent q
#
# example:
#   generalized Hurst exponent for a random vector
#   H=genhurst(np.random.rand(10000,1),3)
#
####################################
# for the generalized Hurst exponent method please refer to:
#
#   T. Di Matteo et al. Physica A 324 (2003) 183-188 
#   T. Di Matteo et al. Journal of Banking & Finance 29 (2005) 827-851
#   T. Di Matteo Quantitative Finance, 7 (2007) 21-36
#
####################################
##   written in Matlab : Tomaso Aste, 30/01/2013 ##
##   translated to Python (3.6) : Peter Rupprecht, p.t.r.rupprecht (AT) gmail.com, 25/05/2017 ##

 
import numpy as np
import warnings
 
def genhurst(S,q):

    L=len(S)       
    if L < 100:
        warnings.warn('Data series very short!')
       
    H = np.zeros((len(range(5,20)),1))
    k = 0
    
    for Tmax in range(5,20):
        
        x = np.arange(1,Tmax+1,1)
        mcord = np.zeros((Tmax,1))
        
        for tt in range(1,Tmax+1):
            dV = S[np.arange(tt,L,tt)] - S[np.arange(tt,L,tt)-tt] 
            VV = S[np.arange(tt,L+tt,tt)-tt]
            N = len(dV) + 1
            X = np.arange(1,N+1,dtype=np.float64)
            Y = VV
            mx = np.sum(X)/N
            SSxx = np.sum(X**2) - N*mx**2
            my = np.sum(Y)/N
            SSxy = np.sum( np.multiply(X,Y))  - N*mx*my
            cc1 = SSxy/SSxx
            cc2 = my - cc1*mx
            ddVd = dV - cc1
            VVVd = VV - np.multiply(cc1,np.arange(1,N+1,dtype=np.float64)) - cc2
            mcord[tt-1] = np.mean( np.abs(ddVd)**q )/np.mean( np.abs(VVVd)**q )
            
        mx = np.mean(np.log10(x))
        SSxx = np.sum( np.log10(x)**2) - Tmax*mx**2
        my = np.mean(np.log10(mcord))
        SSxy = np.sum( np.multiply(np.log10(x),np.transpose(np.log10(mcord)))) - Tmax*mx*my
        H[k] = SSxy/SSxx
        k = k + 1
        
    mH = np.mean(H)/q
    
    return mH
 
