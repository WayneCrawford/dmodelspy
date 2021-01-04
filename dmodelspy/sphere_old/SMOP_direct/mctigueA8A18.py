# Generated with SMOP  0.41-beta
from smop.libsmop import *
# mctigueA8A18.m

    
@function
def mctigueA8A18(tt=None,nu=None,rho=None,zzn=None,*args,**kwargs):
    varargin = mctigueA8A18.varargin
    nargin = mctigueA8A18.nargin

    # first free surface correction
# displacement functions (A8) and (A18) of McTigue (1987)
#==========================================================================
# USGS Software Disclaimer 
# The software and related documentation were developed by the U.S. 
# Geological Survey (USGS) for use by the USGS in fulfilling its mission. 
# The software can be used, copied, modified, and distributed without any 
# fee or cost. Use of appropriate credit is requested.
    
    # The USGS provides no warranty, expressed or implied, as to the correctness 
# of the furnished software or the suitability for any purpose. The software 
# has been tested, but as with any complex software, there could be undetected 
# errors. Users who find errors are requested to report them to the USGS. 
# The USGS has limited resources to assist non-USGS users; however, we make 
# an attempt to fix reported problems and help whenever possible. 
#==========================================================================
    
    sigma=multiply(dot(0.5,tt),exp(- tt))
# mctigueA8A18.m:19
    tau1=copy(sigma)
# mctigueA8A18.m:20
    A8=multiply(multiply(multiply(dot(0.5,sigma),((1 - dot(2,nu)) - multiply(tt,zzn))),exp(multiply(- tt,zzn))),besselj(1,multiply(tt,rho)))
# mctigueA8A18.m:22
    A18=multiply(multiply(multiply(dot(0.5,tau1),(dot(2,(1 - nu)) - multiply(tt,zzn))),exp(multiply(- tt,zzn))),besselj(1,multiply(tt,rho)))
# mctigueA8A18.m:23
    dur=A8 + A18
# mctigueA8A18.m:24