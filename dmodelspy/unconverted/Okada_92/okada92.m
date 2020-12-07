 function [u v w dwdx dwdy eea gamma1 gamma2] = okada92(fault,xi,yi,xf,yf,zt,zb,U,delta,mu,nu,x,y,z)
% Compute the internal displacement, strain and tilt  due to a rectangular 
% dislocation. Based on Okada (1992). All parameters are in SI units 
%
% OUTPUT
% u         horizontal (East component) deformation
% v         horizontal (North component) deformation
% w         vertical (Up component) deformation
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


% [1] Set the parameters for the Okada (1992) dislocation model ***********
% Lame's first parameter
lambda = 2*mu*nu/(1-2*nu);

% change fault geometry units from m to km
xi = 0.001*xi; yi = 0.001*yi; xf = 0.001*xf; yf = 0.001*yf; 
zt = 0.001*zt; zb = 0.001*zb; 

% change bechmarks location and depth units from m to km
x = 0.001*x; y = 0.001*y; z = 0.001*z;

% change dip and strike angles in radians
delta = pi*delta/180; 
phi = atan2((xf-xi),(yf-yi));

% check that dip angle delta is different from zero
if abs(delta) < 1E-8
    error('dip angle delta equal to zero');
end;

% length (L) in KM and width (W) in KM of the rectangular fault
L = sqrt((xf-xi)^2+(yi-yf)^2);
if zb > zt
    W = abs((zb-zt)/sin(delta));
else
    error('fault top zt deeper than fault bottom zb');
end;

% coordinates in KM of left lower corner of the rectangular fault (see Figure 4; Okada, 1992)
x0 = xi + (zb-zt)*cot(delta)*cos(phi);
y0 = yi - (zb-zt)*cot(delta)*sin(phi);

% depth in KM of of left lower corner of the rectangular fault (see Figure 4; Okada, 1992)
c = zb;

% Translate a rotate coordinates in the reference system of Fig. 3 (Okada, 1992)
% Center the coordinates systen in (x0,y0)
xxn = x - x0;
yyn = y - y0; 
% rotate the coordinate system 
xxp  = sin(phi)*xxn  + cos(phi)*yyn;
yyp  = -cos(phi)*xxn + sin(phi)*yyn;
% *************************************************************************


% [2] Compute the displacement and displacement gradient matrix in the ****
%     fault coordinate system 
% first check that the right fault is called
if strcmpi(fault,'strike')==0 && strcmpi(fault,'dip')==0 && strcmpi(fault,'tensile')==0 
    error('Wrong fault name');
end;

% change the sign to the slip U for a strike slip fault so that U >0 for a
% right lateral strike slip fault (same convention of Coulomb 3.1)
if strcmpi(fault,'strike')==1, U = -U; end;

% uA (infinite medium), e.g. Okada, Table 6 and pg 1031
% Chinnery's notation, Okada (1992), equation (15)
  [uA1a uA2a uA3a jA1a jA2a jA3a kA1a kA2a kA3a] = ok92fa(fault,xxp,0,c,mu,lambda,U,delta,yyp,z);
  [uA1b uA2b uA3b jA1b jA2b jA3b kA1b kA2b kA3b] = ok92fa(fault,xxp,W,c,mu,lambda,U,delta,yyp,z);
  [uA1c uA2c uA3c jA1c jA2c jA3c kA1c kA2c kA3c] = ok92fa(fault,xxp-L,0,c,mu,lambda,U,delta,yyp,z);
  [uA1d uA2d uA3d jA1d jA2d jA3d kA1d kA2d kA3d] = ok92fa(fault,xxp-L,W,c,mu,lambda,U,delta,yyp,z);
% uAtilde (infinite medium), e.g. Okada, Table 6 and pg 1031
% Chinnery's notation, Okada (1992), equation (15)
  [uA1ta uA2ta uA3ta jA1ta jA2ta jA3ta kA1ta kA2ta kA3ta] = ok92fa(fault,xxp,0,c,mu,lambda,U,delta,yyp,-z);
  [uA1tb uA2tb uA3tb jA1tb jA2tb jA3tb kA1tb kA2tb kA3tb] = ok92fa(fault,xxp,W,c,mu,lambda,U,delta,yyp,-z);
  [uA1tc uA2tc uA3tc jA1tc jA2tc jA3tc kA1tc kA2tc kA3tc] = ok92fa(fault,xxp-L,0,c,mu,lambda,U,delta,yyp,-z);
  [uA1td uA2td uA3td jA1td jA2td jA3td kA1td kA2td kA3td] = ok92fa(fault,xxp-L,W,c,mu,lambda,U,delta,yyp,-z);
% uB (free surface deformation), e.g. Okada (1992), table 6 and equation (18)
% Chinnery's notation, Okada (1992), equation (15)
  [uB1a uB2a uB3a jB1a jB2a jB3a kB1a kB2a kB3a] = ok92fb(fault,xxp,0,c,mu,lambda,U,delta,yyp,z);
  [uB1b uB2b uB3b jB1b jB2b jB3b kB1b kB2b kB3b] = ok92fb(fault,xxp,W,c,mu,lambda,U,delta,yyp,z);
  [uB1c uB2c uB3c jB1c jB2c jB3c kB1c kB2c kB3c] = ok92fb(fault,xxp-L,0,c,mu,lambda,U,delta,yyp,z);
  [uB1d uB2d uB3d jB1d jB2d jB3d kB1d kB2d kB3d] = ok92fb(fault,xxp-L,W,c,mu,lambda,U,delta,yyp,z);
% uC (depth multiplied term), e.g. Okada (1992), Table 6 and pg 1031
% Chinnery's notation, Okada (1992), equation (15)
  [uC1a uC2a uC3a jC1a jC2a jC3a kC1a kC2a kC3a] = ok92fc(fault,xxp,0,c,mu,lambda,U,delta,yyp,z);
  [uC1b uC2b uC3b jC1b jC2b jC3b kC1b kC2b kC3b] = ok92fc(fault,xxp,W,c,mu,lambda,U,delta,yyp,z);
  [uC1c uC2c uC3c jC1c jC2c jC3c kC1c kC2c kC3c] = ok92fc(fault,xxp-L,0,c,mu,lambda,U,delta,yyp,z);
  [uC1d uC2d uC3d jC1d jC2d jC3d kC1d kC2d kC3d] = ok92fc(fault,xxp-L,W,c,mu,lambda,U,delta,yyp,z);
% *************************************************************************
            

% [3] Compute the displacement (Okada, 1992; Table 6) *********************            
% uA (infinite medium), see also pg 1031
    uA1 = uA1a-uA1b-uA1c+uA1d;
    uA2 = uA2a-uA2b-uA2c+uA2d;
    uA3 = uA3a-uA3b-uA3c+uA3d;
% uAtilde (infinite medium), see also pg 1031
    uA1tilde = uA1ta-uA1tb-uA1tc+uA1td;
    uA2tilde = uA2ta-uA2tb-uA2tc+uA2td;
    uA3tilde = uA3ta-uA3tb-uA3tc+uA3td;
% uB (free surface deformation), see also equation (18)
    uB1 = uB1a-uB1b-uB1c+uB1d;
    uB2 = uB2a-uB2b-uB2c+uB2d;
    uB3 = uB3a-uB3b-uB3c+uB3d;
% uC (depth multiplied term), see also and pg 1031
    uC1 = uC1a-uC1b-uC1c+uC1d;
    uC2 = uC2a-uC2b-uC2c+uC2d;
    uC3 = uC3a-uC3b-uC3c+uC3d;
    
% displacement
Upx =  uA1-uA1tilde + uB1 + z.*uC1;
Upy = (uA2-uA2tilde + uB2 + z.*uC2)*cos(delta) - (uA3-uA3tilde + uB3 + z.*uC3)*sin(delta);
Upz = (uA2-uA2tilde + uB2 - z.*uC2)*sin(delta) + (uA3-uA3tilde + uB3 - z.*uC3)*cos(delta);

% Rotate the horizontal displacement components Upx and Upy back
u = sin(phi)*Upx - cos(phi)*Upy;
v = sin(phi)*Upy + cos(phi)*Upx;
w = Upz;
% *************************************************************************


% [4] Compute the X-derivative (Okada, 1992; Table 7) *********************            
% jA (infinite medium)
    jA1 = jA1a-jA1b-jA1c+jA1d;
    jA2 = jA2a-jA2b-jA2c+jA2d;
    jA3 = jA3a-jA3b-jA3c+jA3d;
% jAtilde (infinite medium)
    jA1tilde = jA1ta-jA1tb-jA1tc+jA1td;
    jA2tilde = jA2ta-jA2tb-jA2tc+jA2td;
    jA3tilde = jA3ta-jA3tb-jA3tc+jA3td;
% jB (free surface deformation) 
    jB1 = jB1a-jB1b-jB1c+jB1d;
    jB2 = jB2a-jB2b-jB2c+jB2d;
    jB3 = jB3a-jB3b-jB3c+jB3d;
% uC (depth multiplied term)
    jC1 = jC1a-jC1b-jC1c+jC1d;
    jC2 = jC2a-jC2b-jC2c+jC2d;
    jC3 = jC3a-jC3b-jC3c+jC3d;
    
% X-derivative
Dpuxx =  jA1-jA1tilde + jB1 + z.*jC1;
Dpuyx = (jA2-jA2tilde + jB2 + z.*jC2)*cos(delta) - (jA3-jA3tilde + jB3 + z.*jC3)*sin(delta);
Dpuzx = (jA2-jA2tilde + jB2 - z.*jC2)*sin(delta) + (jA3-jA3tilde + jB3 - z.*jC3)*cos(delta);
% *************************************************************************


% [5] Compute the Y-derivative (Okada, 1992; Table 8) *********************            
% kA (infinite medium)
    kA1 = kA1a-kA1b-kA1c+kA1d;
    kA2 = kA2a-kA2b-kA2c+kA2d;
    kA3 = kA3a-kA3b-kA3c+kA3d;
% kAtilde (infinite medium)
    kA1tilde = kA1ta-kA1tb-kA1tc+kA1td;
    kA2tilde = kA2ta-kA2tb-kA2tc+kA2td;
    kA3tilde = kA3ta-kA3tb-kA3tc+kA3td;
% kB (free surface deformation) 
    kB1 = kB1a-kB1b-kB1c+kB1d;
    kB2 = kB2a-kB2b-kB2c+kB2d;
    kB3 = kB3a-kB3b-kB3c+kB3d;
% uC (depth multiplied term)
    kC1 = kC1a-kC1b-kC1c+kC1d;
    kC2 = kC2a-kC2b-kC2c+kC2d;
    kC3 = kC3a-kC3b-kC3c+kC3d;
    
% Y-derivative
Dpuxy =  kA1-kA1tilde + kB1 + z.*kC1;
Dpuyy = (kA2-kA2tilde + kB2 + z.*kC2)*cos(delta) - (kA3-kA3tilde + kB3 + z.*kC3)*sin(delta);
Dpuzy = (kA2-kA2tilde + kB2 - z.*kC2)*sin(delta) + (kA3-kA3tilde + kB3 - z.*kC3)*cos(delta);
% *************************************************************************


% [6] Rotate the displacement gradient matrix back ************************
duxx = Dpuxx*sin(phi)^2 - (Dpuxy+Dpuyx)*sin(phi)*cos(phi) + Dpuyy*cos(phi)^2;
duxy = Dpuxy*sin(phi)^2 + (Dpuxx-Dpuyy)*sin(phi)*cos(phi) - Dpuyx*cos(phi)^2;
duyx = Dpuyx*sin(phi)^2 + (Dpuxx-Dpuyy)*sin(phi)*cos(phi) - Dpuxy*cos(phi)^2;
duyy = Dpuyy*sin(phi)^2 + (Dpuxy+Dpuyx)*sin(phi)*cos(phi) + Dpuxx*cos(phi)^2;
dwdx = Dpuzx*sin(phi) - Dpuzy*cos(phi);                                     % ground tilt (East)
dwdy = Dpuzy*sin(phi) + Dpuzx*cos(phi);                                     % ground tilt (North)
% *************************************************************************

% [7] Scale ground tilt ***************************************************
dwdx = 0.001*dwdx;
dwdy = 0.001*dwdy;
% *************************************************************************

% [8] Compute and scale shear and areal strain ****************************
% Strains
eea = 0.001*(duxx + duyy);                                                  % areal strain, equation (5)
gamma1 = 0.001*(duxx - duyy);                                               % shear strain
gamma2 = 0.001*(duxy + duyx);                                               % shear strain
% *************************************************************************