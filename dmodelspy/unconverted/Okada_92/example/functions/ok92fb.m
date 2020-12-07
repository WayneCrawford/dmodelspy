function [ub1 ub2 ub3 jb1 jb2 jb3 kb1 kb2 kb3] = ok92fb(fault,csi,W,c,mu,lambda,U,delta,yyp,z)
% function [ub1 ub2 ub3 jb1 jb2 jb3 kb1 kb2 kb3] = ok92fb(fault,csi,W,c,mu,lambda,U,delta,yyp,z)
% Forward model of fault slip (ub), X-derivative (jb) and Y-derivative (kb) 
% based on Table 6, table 7 and Table 8, column fB, of Okada (1992). 
%
% INPUT
% fault     a string that define the kind of fault: strike, dip or tensile
% csi       dummy coordinate, defined in (12) and used in the integral (13)
% W         width in km of the rectangular fault, see okada92.m
% c         depth in km of of left lower corner of the rectangular fault
% mu        shear modulus
% lambda    first Lame's elastic parameter
% U         fault slip
% delta     dip angle from horizontal reference surface (0.5*pi = vertical fault)
% yyp       rotated coordinate, see okada92.m
% z         depth of internal deformation (km)
%
% OUTPUT
% The variables ub1 ... ka3 corrrespond to the top, middle and bottom equations
% in the second column of Table 6, table 7 and Table 8.
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


% *** GENERAL PARAMETERS **************************************************
eps = 1E-8;  % criteria used to assume that a variable is equal to zero 

% Okada (1992), Table 6 (just below title)
alpha = (lambda+mu)/(lambda+2*mu);
alpha1 = (1-alpha)/alpha;
d = c-z;
p = yyp*cos(delta) + d*sin(delta);
q = yyp*sin(delta) - d*cos(delta);
eta = p-W;                                      % integration variable
R2 = csi.^2 + eta.^2 + q.^2; R = sqrt(R2);                                  
ytilde = eta*cos(delta) + q*sin(delta);                                     
dtilde = eta*sin(delta) - q*cos(delta);                                     
ctilde = dtilde + z;                                                        

% Okada (1992), equation (14)
X11 = 1./(R.*(R+csi)); X32 = (2*R+csi)./(R.^3.*(R+csi).^2); X53 = (8*R.^2+9*R.*csi+3*csi.^2)./(R.^5.*(R+csi).^3);           % csi = xxp
Y11 = 1./(R.*(R+eta)); Y32 = (2*R+eta)./(R.^3.*(R+eta).^2); Y53 = (8*R.^2+9*R.*eta+3*eta.^2)./(R.^5.*(R+eta).^3);
h = q*cos(delta) - z;  Z32 = sin(delta)./R.^3 - h.*Y32;     Z53 = 3*sin(delta)./R.^5 - h.*Y53;
Y0 = Y11 - csi.^2.*Y32; Z0 = Z32 - csi.^2.*Z53;

% THETA (footnote of Table 6)
F0 = zeros(size(q));
c1 = abs(q) <= eps;
c2 = abs(q) > eps;
THETA = c1.*F0 + c2.*atan(csi.*eta./(q.*R));

% check singularity condition (iii), pg 1034, Okada (1992)
F0 = zeros(size(R+csi));
c1 = abs(R+csi) <= eps; 
c2 = abs(R+csi) > eps;
Rcsi = c1.*F0 + c2.*(1./(R+csi));
lnRcsi = -c1.*log(R-csi) + c2.*log(R+csi);

% check singularity condition (iv), pg 1034, Okada (1992)
F0 = zeros(size(R+eta));
c1 = abs(R+eta) <= eps; 
c2 = abs(R+eta) > eps;
Reta = c1.*F0 + c2.*(1./(R+eta));
lnReta = -c1.*log(R-eta) + c2.*log(R+eta);

% I parameters (footnote of Table 6)
F0 = zeros(size(csi));
c1 = abs(csi) <= eps; 
c2 = abs(csi) > eps;
X = sqrt(csi.^2+q.^2);

if cos(delta) < eps                                                         % check singularity for cos(delta)=0
    I3 = 0.5*(eta./(R+dtilde) + ytilde.*q./(R+dtilde).^2 -lnReta); 
    I4 = c1.*F0 + c2.*(0.5*csi.*ytilde./(R+dtilde).^2);                     % check singularity condition (ii), pg 1034, Okada (1992)       
else                          
    I3 = (1/cos(delta))*ytilde./(R+dtilde) - (1/cos(delta)^2)*(lnReta-sin(delta)*log(R+dtilde));
    I4 = c1.*F0 + c2.*((sin(delta)/cos(delta))*csi./(R+dtilde) + (2/cos(delta)^2)*...
                atan2((eta.*(X+q.*cos(delta))+X.*(R+X)*sin(delta)),(csi.*(R+X)*cos(delta))));
end;
    I1 = -(csi./(R+dtilde))*cos(delta) - I4*sin(delta);
    I2 = log(R+dtilde) + I3*sin(delta);

% J and K parameters (footnote of Table 7)
D11 = 1./(R.*(R+dtilde));
    J2 = csi.*ytilde.*D11./(R+dtilde);
    J5 = -(dtilde + ytilde.^2./(R+dtilde)).*D11;
if cos(delta) < eps
    K1 = csi.*q.*D11./(R+dtilde);
    K3 = sin(delta)*(csi.^2.*D11-1)./(R+dtilde);
    J3 = -csi.*(q.^2.*D11-0.5)./(R+dtilde).^2;
    J6 = -ytilde.*(csi.^2.*D11-0.5)./(R+dtilde).^2;
else
    K1 = csi.*(D11-Y11*sin(delta))/cos(delta);
    K3 = (q.*Y11 - ytilde.*D11)/cos(delta);
    J3 = (K1 - J2*sin(delta))/cos(delta);
    J6 = (K3 - J5*sin(delta))/cos(delta);
end;
    J1 = J5*cos(delta) - J6*sin(delta);
    J4 = -csi.*Y11-J2*cos(delta)+J3*sin(delta); 
    K2 = 1./R + K3*sin(delta);
    K4 = csi.*Y11*cos(delta) - K1*sin(delta);

% E-Q parameters (footnote Table 8)
E = sin(delta)./R - ytilde.*q./R.^3;
F = dtilde./R.^3 + csi.^2.*Y32*sin(delta);
G = 2*X11*sin(delta) - ytilde.*q.*X32;
H = dtilde.*q.*X32 + csi.*q.*Y32*sin(delta);
P = cos(delta)./R.^3 + q.*Y32*sin(delta);
Q = 3*ctilde.*dtilde./R.^5 - (z.*Y32+Z32+Z0)*sin(delta);
% *************************************************************************


if strcmpi(fault,'strike')==1
    % STRIKE SLIP displacement: Okada (1992), Table 6, Type Strike
    ub1 = -csi.*q.*Y11 -THETA -alpha1.*I1*sin(delta);                           % first row, second column
    ub2 = -q./R +alpha1.*ytilde*sin(delta)./(R+dtilde);                         % second row, second column
    ub3 = q.^2.*Y11 -alpha1.*I2*sin(delta);                                     % third row, second column
% *************************************************************************
    % STRIKE SLIP X-derivative: Okada (1992), Table 7, Type Tensile
    jb1 =  csi.^2.*q.*Y32 -alpha1*J1*sin(delta);                                % first row, first column
    jb2 =  csi.*q./R.^3   -alpha1*J2*sin(delta);                                  % second row, first column
    jb3 = -csi.*q.^2.*Y32 -alpha1*J3*sin(delta);                                % third row, first column
% *************************************************************************
    % STRIKE SLIP Y-Derivative: Okada (1992), Table 8, Type Tensile
    kb1 = -csi.*F -dtilde.*X11 +alpha1*(csi.*Y11+J4)*sin(delta);                % first row, first column
    kb2 = -E                   +alpha1*(1./R    +J5)*sin(delta);                % second row, first column
    kb3 = q.*F                 -alpha1*(q.*Y11  -J6)*sin(delta);                % third row, first column
% *************************************************************************
elseif strcmpi(fault,'dip')==1
    % DIP SLIP displacement: Okada (1992), Table 6, Type Dip
    ub1 = -q./R +((1-alpha)/alpha).*I3*sin(delta)*cos(delta);                               % first row, second column
    ub2 = -eta.*q.*X11-THETA -((1-alpha)/alpha).*(csi./(R+dtilde))*sin(delta)*cos(delta);   % second row, second column
    ub3 = q.^2.*X11 +((1-alpha)/alpha).*I4*sin(delta)*cos(delta);                           % third row, second column
% *************************************************************************
    % DIP SLIP X-Derivative: Okada (1992), Table 7, Type Dip
    jb1 = csi.*q./R.^3         +alpha1*J4*sin(delta)*cos(delta);                       % first row, first column
    jb2 = eta.*q./R.^3 +q.*Y11 +alpha1*J5*sin(delta)*cos(delta);              % second row, first column
    jb3 = -q.^2./R.^3          +alpha1*J6*sin(delta)*cos(delta);                        % third row, first column
% *************************************************************************
    % DIP SLIP Y-Derivative: Okada (1992), Table 8, Type Dip
    kb1 = -E                           +alpha1*J1*sin(delta)*cos(delta);             % first row, first column
    kb2 = -eta.*G -csi.*Y11*sin(delta) +alpha1*J2*sin(delta)*cos(delta);             % second row, first column
    kb3 = q.*G                         +alpha1*J3*sin(delta)*cos(delta);             % third row, first column
% *************************************************************************
else
    % TENSILE displacement: Okada (1992), Table 6, Type Tensile
    ub1 = q.^2.*Y11 -alpha1.*I3*sin(delta).^2;                                  % first row, second column
    ub2 = q.^2.*X11 +alpha1.*(csi./(R+dtilde))*sin(delta).^2;                   % second row, second column
    ub3 = q.*(eta.*X11+csi.*Y11)-THETA-alpha1.*I4*sin(delta).^2;                % third row, second column
% *************************************************************************
    % TENSILE X-derivative: Okada (1992), Table 7, Type Tensile
    jb1 = -csi.*q.^2.*Y32 -alpha1*J4*sin(delta).^2;                            % first row, first column
    jb2 = -q.^2./R.^3     -alpha1*J5*sin(delta).^2;                                % second row, first column
    jb3 =  q.^3.*Y32      -alpha1*J6*sin(delta).^2;                                  % third row, first column
% *************************************************************************
    % TENSILE Y-derivative: Okada (1992), Table 8, Type Tensile
    kb1 =  q.*F -alpha1*J1*sin(delta).^2;                                       % first row, first column
    kb2 =  q.*G -alpha1*J2*sin(delta).^2;                                       % second row, first column
    kb3 = -q.*H -alpha1*J3*sin(delta).^2;                                       % third row, first column
% *************************************************************************
end


% *************************************************************************
% Displacement: Okada (1992), Table 6
ub1 = (0.5*U/pi)*ub1;
ub2 = (0.5*U/pi)*ub2;
ub3 = (0.5*U/pi)*ub3;
% *************************************************************************
% X-derivative: Okada (1992), Table 7
jb1 = (0.5*U/pi)*jb1;
jb2 = (0.5*U/pi)*jb2;
jb3 = (0.5*U/pi)*jb3;
% *************************************************************************
% Y-derivative: Okada (1992), Table 8, Type Tensile
kb1 = (0.5*U/pi)*kb1;
kb2 = (0.5*U/pi)*kb2;
kb3 = (0.5*U/pi)*kb3;
% *************************************************************************