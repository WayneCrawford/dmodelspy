# Generated with SMOP  0.41-beta
from smop.libsmop import *
# mctigueA7A17.m

    
@function
def mctigueA7A17(tt=None,nu=None,rho=None,zeta=None,*args,**kwargs):
    varargin = mctigueA7A17.varargin
    nargin = mctigueA7A17.nargin

    # first free surface correction
# displacement functions (A7) and (A17) of McTigue (1987)
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
    
    R=sqrt(rho ** 2 + zeta ** 2)
# mctigueA7A17.m:19
    sigma=multiply(dot(0.5,tt),exp(- tt))
# mctigueA7A17.m:20
    tau1=copy(sigma)
# mctigueA7A17.m:21
    A7=multiply(multiply(multiply(dot(0.5,sigma),(dot(2,(1 - nu)) - multiply(tt,zeta))),exp(multiply(tt,zeta))),besselj(0,multiply(tt,R)))
# mctigueA7A17.m:23
    A17=multiply(multiply(multiply(dot(0.5,tau1),((1 - dot(2,nu)) - multiply(tt,zeta))),exp(multiply(tt,zeta))),besselj(0,multiply(tt,R)))
# mctigueA7A17.m:24
    duz=A7 + A17
# mctigueA7A17.m:25