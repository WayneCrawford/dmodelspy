function [u v w dwdx dwdy eea gamma1 gamma2] = okada85fault(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x,y)
% This script approximates a rectangular fault by a matrix of rectangular patches.
% The model uses the rectangular dislocation Green's function by Okada (1985), 
% implemented in the script okada85.m   
%
% All variable/parameters are in SI units, angles are in degrees
% 
% OUTPUT
% u         horizontal (East component) deformation     [m]
% v         horizontal (North component) deformation    [m]
% w         vertical (Up component) deformation         [m]
% dwdx      ground tilt (East component)                [-]
% dwdy      ground tilt (North component)               [-]
% eea       areal strain (du/dx + dv/dy)                [-]
% gamma1    shear strain (du/dx - dv/dy)                [-]
% gamma2    shear strain (du/dy + dv/dx)                [-]
%
% FAULT PARAMETER
% fault     a string that define the kind of fault: strike, dip or tensile
%
% PATCH PARAMETERS (all column vectors with the same length) 
% xi        x start                                     [m]
% yi        y start                                     [m]
% xf        x finish                                    [m]
% yf        y finish                                    [m]
% zt        top                                         [m]
%           (positive downward and defined as depth below the reference surface)
% zb        bottom; zb > zt                             [m]
%           (positive downward and defined as depth below the reference surface)
% U         fault slip                                  [m]
%           strike slip fault: U > 0 right lateral strike slip
%           dip slip fault   : U > 0 reverse slip
%           tensile fault    : U > 0 tensile opening fault
% delta     dip angle from horizontal reference surface (90° = vertical fault)
%           delta can be between 0° and 90° but must be different from zero!
% phi       strike angle from North (90° = fault trace parallel to the East)
%
% CRUST PARAMETERS                                                                  
% mu        shear modulus                               [Pa]
% nu        Poisson's ratio                             [-]
%
% GEODETIC BENCHMARKS
% x,y       benchmark location (must be COLUMN vectors) [m]
% *************************************************************************
% Y. Okada (1985). Surface deformation due to shear and tensile faults in a
% half-space. Bull Seism Soc Am 75(4), 1135-1154.
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

N = length(xi);                                                             % number of patches

for i=1:N                                                                   % loop over patches
    [Up Vp Wp Dwdxp Dwdyp Eeap Gamma1p Gamma2p] = okada85(fault,xi(i),yi(i),xf(i),yf(i),zt(i),zb(i),U(i),delta(i),mu,nu,x,y);
    up(:,:,i) = Up;                                                         % patch contribution to horizontal (East component) deformation
    vp(:,:,i) = Vp;                                                         % patch contribution to horizontal (North component) deformation
    wp(:,:,i) = Wp;                                                         % patch contribution to vertical (Up component) deformation
    dwdxp(:,:,i) = Dwdxp;                                                   % patch contribution to ground tilt (East component)
    dwdyp(:,:,i) = Dwdyp;                                                   % patch contribution to ground tilt (North component)
    eeap(:,:,i) = Eeap;                                                     % patch contribution to areal strain
    gamma1p(:,:,i) = Gamma1p;                                               % patch contribution to shear strain
    gamma2p(:,:,i) = Gamma2p;                                               % patch contribution to ahear strain
end;

u = sum(up,3);                                                              % horizontal (East component) deformation
v = sum(vp,3);                                                              % horizontal (North component) deformation
w = sum(wp,3);                                                              % vertical (Up component) deformation
dwdx = sum(dwdxp,3);                                                        % ground tilt (East component)
dwdy = sum(dwdyp,3);                                                        % ground tilt (North component)
eea = sum(eeap,3);                                                          % areal strain
gamma1 = sum(gamma1p,3);                                                    % shear strain
gamma2 = sum(gamma2p,3);                                                    % shear strain
