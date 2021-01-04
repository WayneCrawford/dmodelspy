function [eea gamma1 gamma2] = mctigue3Dstrain(x0,y0,z0,P_G,a,nu,x,y,z)
% 3D Green's function for spherical source
% 
% all parameters are in SI (MKS) units
% eea       areal strain
% gamma1    shear strain
% gamma2    shear strain
% x0,y0     coordinates of the center of the sphere 
% z0        depth of the center of the sphere (positive downward and
%              defined as distance below the reference surface)
% P_G       dimensionless excess pressure (pressure/shear modulus)
% a         radius of the sphere
% nu        Poisson's ratio
% x,y       benchmark location
% z         depth within the crust (z=0 is the free surface)
% 
% Reference ***************************************************************
% Battaglia M. et al. Implementing the spherical source
% equation (5) and (6)
% *************************************************************************
% USGS Software Disclaimer 
% The software and related documentation were developed by the U.S. 
% Geological Survey (USGS) for use by the USGS in fulfilling its mission. 
% The software can be used, copied, modified, and distributed without any 
% fee or cost. Use of appropriate credit is requested. 
%
% The USGS provides no warranty, expressed or implied, as to the correctness 
% of the furnished software or the suitability for any purpose. The software 
% has been tested, but as with any complex software, there could be undetected 
% errors. Users who find errors are requested to report them to the USGS. 
% The USGS has limited resources to assist non-USGS users; however, we make 
% an attempt to fix reported problems and help whenever possible. 
%==========================================================================


% STRAIN ******************************************************************
h = 0.001*abs(max(x)-min(x));                                              % finite difference step

% derivatives
[up , ~, ~] = mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x+h,y,z);
[um , ~, ~] = mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x-h,y,z);
dudx = 0.5*(up - um)/h;

[up , ~, ~] = mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x,y+h,z);
[um , ~, ~] = mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x,y-h,z);
dudy = 0.5*(up - um)/h;

[~, vp , ~] = mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x+h,y,z);
[~, vm , ~] = mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x-h,y,z);
dvdx = 0.5*(vp - vm)/h;

[~, vp , ~] = mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x,y+h,z);
[~, vm , ~] = mctigue3Ddispl(x0,y0,z0,P_G,a,nu,x,y-h,z);
dvdy = 0.5*(vp - vm)/h;

% Strains
eea = dudx + dvdy;                                                          % areal strain, equation (5)
gamma1 = dudx - dvdy;                                                       % shear strain
gamma2 = dudy + dvdx;                                                       % shear strain





