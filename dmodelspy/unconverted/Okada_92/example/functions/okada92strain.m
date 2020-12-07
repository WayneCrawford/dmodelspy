function [eea gamma1 gamma2] = okada92strain(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x,y,z)
% function [u v w] = okada92strain(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x,y,z)
% Rectangular dislocation forward model
% Compute the internal strain due to a rectangular dislocation.
% Based on Okada (1992). All parameters are in SI units 
%
% OUTPUT
% eea       areal strain
% gamma1    shear strain
% gamma2    shear strain
%
% FAULT PARAMETERS
% fault     a string that define the kind of fault: strike, dip or tensile
% xi        x start
% yi        y start
% xf        x finish
% yf        y finish
% zt        top
%           (positive downward and defined as depth below the reference surface)
% zb        bottom; zb > zt
%           (positive downward and defined as depth below the reference surface)
% U         fault slip
%           strike slip fault: U > 0 right lateral strike slip
%           dip slip fault   : U > 0 reverse slip
%           tensile fault    : U > 0 tensile opening fault
% delta     dip angle from horizontal reference surface (90° = vertical fault)
%           delta can be between 0° and 90° but must be different from zero!
% phi       strike angle from North (90° = fault trace parallel to the East)
%
% CRUST PARAMETERS                                                                  
% mu        shear modulus
% nu        Poisson's ratio
%
% GEODETIC BENCHMARKS
% x,y       benchmark location (must be COLUMN vectors)
% z         depth of internal deformation (z=0 is the free surface)
% *************************************************************************
% Y. Okada (1992). Internal deformation due to shear and tensile faults in 
% a half-space. Bull Seism Soc Am 82(2), 1018-1040.
% *************************************************************************
%==========================================================================
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


h = 1E-6*abs(max(x)-min(x));                                              % finite difference step

% derivatives
[UP , ~, ~] = okada92fault(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x+h,y,z); 
[UM , ~, ~] = okada92fault(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x-h,y,z);
dudx = 0.5*(UP - UM)/h;

[up , ~, ~] = okada92fault(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x,y+h,z);
[um , ~, ~] = okada92fault(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x,y-h,z);
dudy = 0.5*(up - um)/h;

[~, vp , ~] = okada92fault(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x+h,y,z); 
[~, vm , ~] = okada92fault(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x-h,y,z);
dvdx = 0.5*(vp - vm)/h;

[~, VP , ~] = okada92fault(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x,y+h,z);
[~, VM , ~] = okada92fault(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x,y-h,z);
dvdy = 0.5*(VP - VM)/h;

% Strains
eea = dudx + dvdy;                                                          % areal strain, equation (5)
gamma1 = dudx - dvdy;                                                       % shear strain
gamma2 = dudy + dvdx;                                                       % shear strain