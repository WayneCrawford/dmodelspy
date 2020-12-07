# Generated with SMOP  0.41-beta
from smop.libsmop import *
# mctigue3Dstrain.m

    
@function
def mctigue3Dstrain(x0=None,y0=None,z0=None,P_G=None,a=None,nu=None,x=None,y=None,z=None,*args,**kwargs):
    varargin = mctigue3Dstrain.varargin
    nargin = mctigue3Dstrain.nargin

    # 3D Green's function for spherical source
# 
# all parameters are in SI (MKS) units
# eea       areal strain
# gamma1    shear strain
# gamma2    shear strain
# x0,y0     coordinates of the center of the sphere 
# z0        depth of the center of the sphere (positive downward and
#              defined as distance below the reference surface)
# P_G       dimensionless excess pressure (pressure/shear modulus)
# a         radius of the sphere
# nu        Poisson's ratio
# x,y       benchmark location
# z         depth within the crust (z=0 is the free surface)
# 
# Reference ***************************************************************
# Battaglia M. et al. Implementing the spherical source
# equation (5) and (6)
# *************************************************************************
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
    
    # STRAIN ******************************************************************
    h=dot(0.001,abs(max(x) - min(x)))
# mctigue3Dstrain.m:37
    
    # derivatives
    up,__,__=mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x + h,y,z,nargout=3)
# mctigue3Dstrain.m:40
    um,__,__=mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x - h,y,z,nargout=3)
# mctigue3Dstrain.m:41
    dudx=dot(0.5,(up - um)) / h
# mctigue3Dstrain.m:42
    up,__,__=mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x,y + h,z,nargout=3)
# mctigue3Dstrain.m:44
    um,__,__=mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x,y - h,z,nargout=3)
# mctigue3Dstrain.m:45
    dudy=dot(0.5,(up - um)) / h
# mctigue3Dstrain.m:46
    __,vp,__=mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x + h,y,z,nargout=3)
# mctigue3Dstrain.m:48
    __,vm,__=mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x - h,y,z,nargout=3)
# mctigue3Dstrain.m:49
    dvdx=dot(0.5,(vp - vm)) / h
# mctigue3Dstrain.m:50
    __,vp,__=mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x,y + h,z,nargout=3)
# mctigue3Dstrain.m:52
    __,vm,__=mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x,y - h,z,nargout=3)
# mctigue3Dstrain.m:53
    dvdy=dot(0.5,(vp - vm)) / h
# mctigue3Dstrain.m:54
    # Strains
    eea=dudx + dvdy
# mctigue3Dstrain.m:57
    
    gamma1=dudx - dvdy
# mctigue3Dstrain.m:58
    
    gamma2=dudy + dvdx
# mctigue3Dstrain.m:59
    