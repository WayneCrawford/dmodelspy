# Generated with SMOP  0.41-beta
from smop.libsmop import *
# mctigue2D_test.m

    # test of the 2D finite spherical source
    
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
    
#clear('all')
#close_('all')
#clc
#scrsz=get(0,'ScreenSize')
# mctigue2D_test.m:18
# read input data for flat half-space model (McTigue, 1987)
__,z0=textread('inputdata.txt','%s %f',1,nargout=2)
# mctigue2D_test.m:21

__,P_G=textread('inputdata.txt','%s %f',1,'headerlines',1,nargout=2)
# mctigue2D_test.m:22

__,a=textread('inputdata.txt','%s %f',1,'headerlines',2,nargout=2)
# mctigue2D_test.m:23

__,nu=textread('inputdata.txt','%s %f',1,'headerlines',3,nargout=2)
# mctigue2D_test.m:24

# load FEM model
Xu,U=textread('FEM/udispl-0m.txt','%f %f','headerlines',1,nargout=2)
# mctigue2D_test.m:27
Xv,V=textread('FEM/vdispl-0m.txt','%f %f','headerlines',1,nargout=2)
# mctigue2D_test.m:28
Xw,W=textread('FEM/wdispl-0m.txt','%f %f','headerlines',1,nargout=2)
# mctigue2D_test.m:29
# DISPLACEMENT ************************************************************
# compute deformation at the surface
z=0
# mctigue2D_test.m:33
x=(arange(1e-07,5000,100))
# mctigue2D_test.m:33
y=copy(x)
# mctigue2D_test.m:33

u3,v3,w3=mctigue3Ddispl(0,0,z0,P_G,a,nu,x,y,z,nargout=3)
# mctigue2D_test.m:34

u,v,w,dwdx,dwdy=mctigue2D(0,0,z0,P_G,a,nu,x,y,nargout=5)
# mctigue2D_test.m:35

# compare 3D, 2D and FEM models
figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
subplot(1,3,1)
plot(x,u3,x,u,'r.',Xu,U,'g--')
xlabel('x (m)')
ylabel('x-displacement (m) at z=0')
subplot(1,3,2)
plot(x,v3,x,v,'r.',Xv,V,'g--')
xlabel('x (m)')
ylabel('y-displacement (m) at z=0')
subplot(1,3,3)
plot(x,w3,x,w,'r.',Xw,W,'g--')
xlabel('x (m)')
ylabel('z-displacement (m) at z=0')
legend('3D','2D','FEM')
legend('boxoff')
#filename = 'compare2D3DFEM'; print('-djpeg',filename);

# GROUND TILT ************************************************************* 
# load FE model
X1,dWdx=textread('FEM/dwdx-0m.txt','%f %f','headerlines',1,nargout=2)
# mctigue2D_test.m:48
X2,dWdy=textread('FEM/dwdy-0m.txt','%f %f','headerlines',1,nargout=2)
# mctigue2D_test.m:49
# compare ground tilt against FE model
figure('Position',concat([10,scrsz(4) / 5,scrsz(3) / 3,scrsz(4) / 4]))
subplot(1,2,1)
plot(x,dwdx,'r.',X1,dWdx,'g--')
xlabel('x (m)')
ylabel('ground tilt (E)')
legend('2D','FEM')
legend('boxoff')
subplot(1,2,2)
plot(x,dwdy,'r.',X2,dWdy,'g--')
xlabel('x (m)')
ylabel('ground tilt (N)')
#filename = 'comparetiltFEM'; print('-djpeg',filename);

# PRINT TABLE *************************************************************
fid=fopen('Table1.txt','w+')
# mctigue2D_test.m:59
M=concat([x.T / 1000,y.T / 1000,dot(1000,u.T),dot(1000,v.T),dot(1000,w.T),dot(1000000.0,dwdx.T),dot(1000000.0,dwdy.T),dot(1000,u3.T),dot(1000,v3.T),dot(1000,w3.T)])
# mctigue2D_test.m:60
fprintf(fid,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',M.T)