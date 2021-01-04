# Generated with SMOP  0.41-beta
from smop.libsmop import *
# mctigueA7A17a.m

    
@function
def mctigueA7A17a(tt=None,nu=None,rho=None,zeta=None,*args,**kwargs):
    varargin = mctigueA7A17a.varargin
    nargin = mctigueA7A17a.nargin

    # sixth order free surface correction
# displacement functions (A7) and (A17) of McTigue (1987)
# Hankel transforms from Tellez et al (1997). Tables of Fourier, Laplace
# and Hankel transforms of n-dimensional generalized function. Acta 
# Applicandae Mathematicae, 48, 235-284 
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
# mctigueA7A17a.m:23
    # equation (48)
    sigma=multiply(dot(1.5,(tt + tt ** 2)),exp(- tt)) / (7 - dot(5,nu))
# mctigueA7A17a.m:26
    # missing Hankel transform of second part of equation (49)
    tau=dot(tt ** 2.0,exp(- tt)) / (7 - dot(5,nu))
# mctigueA7A17a.m:28
    A7=multiply(multiply(multiply(dot(0.5,sigma),(dot(2,(1 - nu)) - multiply(tt,zeta))),exp(multiply(tt,zeta))),besselj(0,multiply(tt,R)))
# mctigueA7A17a.m:30
    A17=multiply(multiply(multiply(dot(0.5,tau),((1 - dot(2,nu)) - multiply(tt,zeta))),exp(multiply(tt,zeta))),besselj(0,multiply(tt,R)))
# mctigueA7A17a.m:31
    duz=A7 + A17
# mctigueA7A17a.m:32