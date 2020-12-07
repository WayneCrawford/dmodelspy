function [u v w dwdx dwdy eea gamma1 gamma2 duxx duxy duyx duyy] = okada85(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x,y)
% rectangular dislocation Green's function (forward model)
% compute the displacement, strain and tilt due to a rectangular
% dislocation based on Okada (1985). All parameters are in SI units 
%
% OUTPUT
% u         horizontal (East component) deformation
% v         horizontal (North component) deformation
% w         vertical (Up component) deformation
% dwdx      ground tilt (East component)
% dwdy      ground tilt (North component)
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
% *************************************************************************
% Y. Okada (1985). Surface deformation due to shear and tensile faults in a
% half-space. Bull Seism Soc Am 75(4), 1135-1154.
% *************************************************************************
% Note: 
% the displacement derivatives duxx duxy duyx duyy are presented only
% for testing purposes (check against Table 2 of Okada, 1985). The
% dimension of duxx duxy duyx duyy are m/km
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


% [1] Set the parameters for the Okada (1985) dislocation model ***********
% Lame's first parameter
lambda = 2*mu*nu/(1-2*nu);

% change fault geometry units from m to km
xi = 0.001*xi; yi = 0.001*yi; xf = 0.001*xf; yf = 0.001*yf; 
zt = 0.001*zt; zb = 0.001*zb; 

% change bechmarks location units from m to km
x = 0.001*x; y = 0.001*y;

% change dip and strike angles in radians
delta = pi*delta/180; 
phi = atan2((xf-xi),(yf-yi));

% check that dip angle delta is different from zero
if abs(delta) < 1E-16
    error('dip angle delta equal to zero');
end;

% L, W      length (L) in KM and width (W) in KM of the rectangular fault
L = sqrt((xf-xi)^2+(yi-yf)^2);
if zb > zt
    W = abs((zb-zt)/sin(delta));
else
    error('fault top zt deeper than fault bottom zb');
end;

% x0,y0     coordinates in KM of left lower corner of the rectangular fault 
%                                               (see Figure 1; Okada, 1985)
x0 = xi + (zb-zt)*cot(delta)*cos(phi);
y0 = yi - (zb-zt)*cot(delta)*sin(phi);

% z0        depth in KM of of left lower corner of the rectangular fault 
z0 = zb;

% translate the coordinates of the points where the displacement is computed
% in the coordinates systen centered in (x0,y0)
xxn = x - x0;
yyn = y - y0; 
d = z0;

% rotate the coordinate system to be coherent with the model coordinate
% system of Figure 1 (Okada, 1985)
xxp  = sin(phi)*xxn  + cos(phi)*yyn;
yyp  = -cos(phi)*xxn + sin(phi)*yyn;
% *************************************************************************

% [2] Compute the displacement and displacement gradient matrix ***********
% Okada (1985), equation (30)
p = yyp*cos(delta) + d*sin(delta);
q = yyp*sin(delta) - d*cos(delta);

% first check that the right fault is called
if strcmpi(fault,'strike')==0 && strcmpi(fault,'dip')==0 && strcmpi(fault,'tensile')==0 
    error('Wrong fault name');
end;

if strcmpi(fault,'strike')==1
    % strike slip fault displacement and displacement gradient matrix. Okada (1985), equation (25)
    [ux1 uy1 uz1 duxx1 duxy1 duyx1 duyy1 duzx1 duzy1] = ok8525(xxp,p,q,mu,lambda,-U,delta);
    [ux2 uy2 uz2 duxx2 duxy2 duyx2 duyy2 duzx2 duzy2] = ok8525(xxp,p-W,q,mu,lambda,-U,delta);
    [ux3 uy3 uz3 duxx3 duxy3 duyx3 duyy3 duzx3 duzy3] = ok8525(xxp-L,p,q,mu,lambda,-U,delta);
    [ux4 uy4 uz4 duxx4 duxy4 duyx4 duyy4 duzx4 duzy4] = ok8525(xxp-L,p-W,q,mu,lambda,-U,delta);
elseif strcmpi(fault,'dip')==1
    % dip slip fault displacement and displacement gradient matrix. Okada (1985), equation (26)
    [ux1 uy1 uz1 duxx1 duxy1 duyx1 duyy1 duzx1 duzy1] = ok8526(xxp,p,q,mu,lambda,U,delta);
    [ux2 uy2 uz2 duxx2 duxy2 duyx2 duyy2 duzx2 duzy2] = ok8526(xxp,p-W,q,mu,lambda,U,delta);
    [ux3 uy3 uz3 duxx3 duxy3 duyx3 duyy3 duzx3 duzy3] = ok8526(xxp-L,p,q,mu,lambda,U,delta);
    [ux4 uy4 uz4 duxx4 duxy4 duyx4 duyy4 duzx4 duzy4] = ok8526(xxp-L,p-W,q,mu,lambda,U,delta);
else 
    % tensile fault displacement and displacement gradient matrix. Okada (1985), equation (27)
    [ux1 uy1 uz1 duxx1 duxy1 duyx1 duyy1 duzx1 duzy1] = ok8527(xxp,p,q,mu,lambda,U,delta);
    [ux2 uy2 uz2 duxx2 duxy2 duyx2 duyy2 duzx2 duzy2] = ok8527(xxp,p-W,q,mu,lambda,U,delta);
    [ux3 uy3 uz3 duxx3 duxy3 duyx3 duyy3 duzx3 duzy3] = ok8527(xxp-L,p,q,mu,lambda,U,delta);
    [ux4 uy4 uz4 duxx4 duxy4 duyx4 duyy4 duzx4 duzy4] = ok8527(xxp-L,p-W,q,mu,lambda,U,delta);    
end;

% displacement, Chinnery's notation, Okada (1985), equation (24)
Upx = (ux1 - ux2 - ux3 + ux4);
Upy = (uy1 - uy2 - uy3 + uy4);

% displacement gradient matrix, Chinnery's notation, Okada (1985), equation (24)
Dpuxx = (duxx1 - duxx2 - duxx3 + duxx4);
Dpuxy = (duxy1 - duxy2 - duxy3 + duxy4);
Dpuyx = (duyx1 - duyx2 - duyx3 + duyx4);
Dpuyy = (duyy1 - duyy2 - duyy3 + duyy4);
Dpuzx = (duzx1 - duzx2 - duzx3 + duzx4);
Dpuzy = (duzy1 - duzy2 - duzy3 + duzy4);

% Rotate the horizontal displacement components Upx and Upy back
u = sin(phi)*Upx - cos(phi)*Upy;
v = sin(phi)*Upy + cos(phi)*Upx;
w = (uz1 - uz2 - uz3 + uz4);

% Rotate the displacement gradient matrix back
duxx = Dpuxx*sin(phi)^2 - (Dpuxy+Dpuyx)*sin(phi)*cos(phi) + Dpuyy*cos(phi)^2;
duxy = Dpuxy*sin(phi)^2 + (Dpuxx-Dpuyy)*sin(phi)*cos(phi) - Dpuyx*cos(phi)^2;
duyx = Dpuyx*sin(phi)^2 + (Dpuxx-Dpuyy)*sin(phi)*cos(phi) - Dpuxy*cos(phi)^2;
duyy = Dpuyy*sin(phi)^2 + (Dpuxy+Dpuyx)*sin(phi)*cos(phi) + Dpuxx*cos(phi)^2;
dwdx = Dpuzx*sin(phi) - Dpuzy*cos(phi);                                     % ground tilt (East)
dwdy = Dpuzy*sin(phi) + Dpuzx*cos(phi);                                     % ground tilt (North)
% *************************************************************************


% [3] Scale ground tilt ***************************************************
dwdx = 0.001*dwdx;
dwdy = 0.001*dwdy;

% [4] Compute and scale shear and areal strain ****************************
% Strains
eea = 0.001*(duxx + duyy);                                                  % areal strain, equation (5)
gamma1 = 0.001*(duxx - duyy);                                               % shear strain
gamma2 = 0.001*(duxy + duyx);                                               % shear strain
% *************************************************************************