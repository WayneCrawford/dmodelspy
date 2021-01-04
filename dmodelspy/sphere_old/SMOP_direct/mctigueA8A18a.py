# Generated with SMOP  0.41-beta
from smop.libsmop import *
# mctigueA8A18a.m

    
@function
def mctigueA8A18a(tt=None,nu=None,rho=None,zeta=None,*args,**kwargs):
    varargin = mctigueA8A18a.varargin
    nargin = mctigueA8A18a.nargin

    # sixth order free surface correction
# displacement functions (A8) and (A18) of McTigue (1987)
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
    
    # equation (48)
    sigma=dot(tt ** 2.0,exp(- tt)) / (7 - dot(5,nu))
# mctigueA8A18a.m:23
    # missing Hankel transform of second part of equation (49)
    tau=dot(tt ** 2.0,exp(- tt)) / (7 - dot(5,nu))
# mctigueA8A18a.m:25
    A8=multiply(multiply(multiply(dot(0.5,sigma),((1 - dot(2,nu)) - multiply(tt,zeta))),exp(multiply(- tt,zeta))),besselj(1,multiply(tt,rho)))
# mctigueA8A18a.m:28
    A18=multiply(multiply(multiply(dot(0.5,tau),(dot(2,(1 - nu)) - multiply(tt,zeta))),exp(multiply(- tt,zeta))),besselj(1,multiply(tt,rho)))
# mctigueA8A18a.m:29
    dur=A8 + A18
# mctigueA8A18a.m:30