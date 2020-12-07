# Generated with SMOP  0.41-beta
from smop.libsmop import *
# mctigue3Ddispl.m

    
@function
def mctigue3D(x0=None,y0=None,z0=None,P_G=None,a=None,nu=None,x=None,y=None,z=None,*args,**kwargs):
    varargin = mctigue3D.varargin
    nargin = mctigue3D.nargin

    # 3D Green's function for spherical source
# 
# all parameters are in SI (MKS) units
# u         horizontal (East component) deformation
# v         horizontal (North component) deformation
# w         vertical (Up component) deformation
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
# McTigue, D.F. (1987). Elastic Stress and Deformation Near a Finite 
# Spherical Magma Body: Resolution of the Point Source Paradox. J. Geophys.
# Res. 92, 12,931-12,940.
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
    
    # translate the coordinates of the points where the displacement is computed
# in the coordinates system centered in (x0,y0)
    xxn=x - x0
# mctigue3Ddispl.m:39
    yyn=y - y0
# mctigue3Ddispl.m:40
    # radial distance from source center to points where we compute ur and uz
    r=sqrt(xxn ** 2 + yyn ** 2)
# mctigue3Ddispl.m:43
    # dimensionless parameters used in the formulas
    rho=r / z0
# mctigue3Ddispl.m:46
    e=a / z0
# mctigue3Ddispl.m:46
    zeta=z / z0
# mctigue3Ddispl.m:46
    # DIMENSIONLESS DISPLACEMENT **********************************************
    
    # leading order solution for a pressurized sphere in an unbounded region
# based on equation (11) of McTigue (1987)
    uz0=dot(dot(e ** 3,0.25),(1 - zeta)) / (rho ** 2 + (1 - zeta) ** 2) ** 1.5
# mctigue3Ddispl.m:53
    
    #uz0 = e^3*0.25./(rho.^2+(1-zeta).^2).^1.5;                           # return equation (14) when zeta=0
    ur0=dot(dot(e ** 3,0.25),rho) / (rho ** 2 + (1 - zeta) ** 2) ** 1.5
# mctigue3Ddispl.m:55
    
    # first free surface correction, equations (A7), (A8), (A17) and (A18) 
# return equation (22) and (23) when csi = 0
    Auz1=zeros(length(zeta),length(rho))
# mctigue3Ddispl.m:60
    Aur1=zeros(length(zeta),length(rho))
# mctigue3Ddispl.m:61
    for i in arange(1,length(zeta)).reshape(-1):
        for j in arange(1,length(rho)).reshape(-1):
            Auz1[i,j]=quadl(lambda xx=None: mctigueA7A17(xx,nu,rho(j),zeta(i)),0,50)
# mctigue3Ddispl.m:64
            Aur1[i,j]=quadl(lambda xx=None: mctigueA8A18(xx,nu,rho(j),zeta(i)),0,50)
# mctigue3Ddispl.m:65
    
    # higher order cavity correction, equations (38) and (39)
    R=sqrt(rho ** 2 + (1 - zeta) ** 2)
# mctigue3Ddispl.m:71
    sint=rho / R
# mctigue3Ddispl.m:72
    cost=(1 - zeta) / R
# mctigue3Ddispl.m:72
    C3=concat([dot(e,(1 + nu)) / (dot(12,(1 - nu))),dot(dot(5,e ** 3),(2 - nu)) / (dot(24,(7 - dot(5,nu))))])
# mctigue3Ddispl.m:73
    
    D3=concat([dot(- e ** 3,(1 + nu)) / 12,dot(e ** 5,(2 - nu)) / (dot(4,(7 - nu)))])
# mctigue3Ddispl.m:74
    
    P0=1
# mctigue3Ddispl.m:75
    P2=dot(0.5,(dot(3,cost ** 2) - 1))
# mctigue3Ddispl.m:75
    
    dP0=0
# mctigue3Ddispl.m:76
    dP2=dot(3,cost)
# mctigue3Ddispl.m:76
    ur38=dot(dot(- 0.5,P0),D3(1)) / R ** 2 + multiply((dot(C3(2),(5 - dot(4,nu))) - dot(1.5,D3(2)) / R ** 2),P2) / R ** 2
# mctigue3Ddispl.m:77
    
    ut39=dot(- (dot(dot(2,C3(1)),(1 - nu)) - dot(0.5,D3(1)) / R ** 2),dP0) - multiply((dot(C3(2),(1 - dot(2,nu))) + dot(0.5,D3(2)) / R ** 2),dP2) / R ** 2
# mctigue3Ddispl.m:78
    
    ut39=multiply(ut39,sint)
# mctigue3Ddispl.m:80
    Auz3=multiply(ur38,cost) - multiply(ut39,sint)
# mctigue3Ddispl.m:81
    
    Aur3=multiply(ur38,sint) + multiply(ut39,cost)
# mctigue3Ddispl.m:82
    
    # sixth order surface correction, return equation (50) and (51) when zeta=0
    Auz6=zeros(length(zeta),length(rho))
# mctigue3Ddispl.m:85
    Aur6=zeros(length(zeta),length(rho))
# mctigue3Ddispl.m:86
    for i in arange(1,length(zeta)).reshape(-1):
        for j in arange(1,length(rho)).reshape(-1):
            Auz6[i,j]=quadl(lambda xx=None: mctigueA7A17a(xx,nu,rho(j),zeta(i)),0,50)
# mctigue3Ddispl.m:89
            Aur6[i,j]=quadl(lambda xx=None: mctigueA8A18a(xx,nu,rho(j),zeta(i)),0,50)
# mctigue3Ddispl.m:90
    
    # total surface displacement, return equation (52) and (53) when zeta = 0
    if size(uz0,1) == size(Auz1,2):
        Auz1=Auz1.T
# mctigue3Ddispl.m:96
        Aur1=Aur1.T
# mctigue3Ddispl.m:97
    
    if size(uz0,1) == size(Auz6,2):
        Auz6=Auz6.T
# mctigue3Ddispl.m:100
        Aur6=Aur6.T
# mctigue3Ddispl.m:101
    
    uz=uz0 + dot(e ** 3,Auz1) + dot(e ** 3,Auz3) + dot(e ** 6,Auz6)
# mctigue3Ddispl.m:103
    #uz = e^3*Auz6;
    ur=ur0 + dot(e ** 3,Aur1) + dot(e ** 3,Aur3) + dot(e ** 6,Aur6)
# mctigue3Ddispl.m:105
    # displacement components
    u=multiply(ur,xxn) / r
# mctigue3Ddispl.m:108
    
    v=multiply(ur,yyn) / r
# mctigue3Ddispl.m:109
    
    w=copy(uz)
# mctigue3Ddispl.m:110
    
    # *************************************************************************
    
    # DIMENSIONAL DISPLACEMENT ************************************************
    u=dot(dot(u,P_G),z0)
# mctigue3Ddispl.m:114
    v=dot(dot(v,P_G),z0)
# mctigue3Ddispl.m:115
    w=dot(dot(w,P_G),z0)
# mctigue3Ddispl.m:116
    # *************************************************************************
    