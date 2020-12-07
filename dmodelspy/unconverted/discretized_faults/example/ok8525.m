function [ux uy uz duxx duxy duyx duyy duzx duzy] = ok8525(csi,eta,q,mu,lambda,U1,delta)
%function [ux uy uz duxx duxy duyx duyy duzx duzy ] = ok8525(csi,eta,q,mu,lambda,U1,delta)
% strike slip fault, Green's function (forward model) based on eq. (25) of
% Okada (1985). All paremeters are in SI units
% csi       dummy coordinate, defined in (22) and used in the integral (23)
% eta       dummy coordinate, defined in (22) and used in the integral (23) 
% q         projected coordinate, see (12) 
% mu        shear modulus
% lambda    Lame's first elastic parameter
% delta     dip angle from horizontal reference surface (90 degrees = vertical fault)
% U1        strike slip (U1 > 0 : left lateral strike slip)
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

% *** GENERAL PARAMETERS **************************************************
alpha = mu/(lambda+mu);

% Okada (1985), equation (30)
ytilde = eta*cos(delta) + q*sin(delta);
dtilde = eta*sin(delta) - q*cos(delta);
R2 = csi.^2 + ytilde.^2 + dtilde.^2; R = sqrt(R2);
X2 = csi.^2 + q.^2; X = sqrt(X2);

% check singularity condition (iii), pg 1148, Okada (1985)
if abs(R+eta) < 1E-16
    Reta = 0;
    lnReta = -log(R-eta);
else
    Reta = 1./(R+eta);
    lnReta = log(R+eta);
end;

% Okada (1985), equation (36)
Ac = (2*R+csi)./(R.^3.*(R+csi).^2);
An = (2*R+eta)./(R.^3.*(R+eta).^2);

% Okada (1985), equation (40) and (41)
% check singularity for cos(delta)=0
if cos(delta) < 1E-16        % Okada (1985), equation (41)
    K1 = alpha*csi.*q./(R.*(R+dtilde).^2);
    K3 = alpha*(sin(delta)./(R+dtilde)).*(csi.^2./(R.*(R+dtilde)) - 1);
else
    K3 = alpha*(1/cos(delta))*(q.*Reta./R - ytilde./(R.*(R+dtilde)));
    K1 = alpha*(csi./cos(delta)).*(1./(R.*(R+dtilde)) - sin(delta)*Reta./R);
end;    
    K2 = alpha*(-sin(delta)./R + q*cos(delta).*Reta./R) -K3;
% *************************************************************************

% *** DISPLACEMENT ********************************************************
% check singularity for cos(delta)=0
if cos(delta) < 1E-16        % Okada (1985), equation (29)
    I1 = -0.5*alpha*csi.*q./(R+dtilde).^2;
    I3 = 0.5*alpha*(eta./(R+dtilde) + ytilde.*q./(R+dtilde).^2 - lnReta);
    I4 = -alpha*q./(R+dtilde);
    I5 = -alpha*csi*sin(delta)./(R+dtilde);
else                        % Okada (1985), equation (28)   
    I4 = alpha*(1/cos(delta))*(log(R+dtilde)-sin(delta)*lnReta);
    if abs(csi) < 1E-16     % check singularity condition (ii), pp 1148, Okada (1985)
        I5 = 0;
    else
        I5 = alpha*(2/cos(delta))*atan((eta.*(X+q.*cos(delta)) + X.*(R+X)*sin(delta))./(csi.*(R+X)*cos(delta)));
    end;
    I3 = alpha*((1/cos(delta))*ytilde./(R+dtilde) - lnReta) + tan(delta)*I4;
    I1 = alpha*((-1/cos(delta))*(csi./(R+dtilde))) - tan(delta)*I5;
end;
I2 = alpha*(-lnReta) - I3;

% strike slip displacement: Okada (1985), equation (25)
if abs(q) < 1E-16           % check singularity condition (i), pp 1148, Okada (1985)
    ux = -(U1/(2*pi)).*(csi.*q.*Reta./R + I1*sin(delta));
else    
    ux = -(U1/(2*pi)).*(csi.*q.*Reta./R + atan(csi.*eta./(q.*R)) + I1*sin(delta));
end;
uy = -(U1/(2*pi)).*(ytilde.*q.*Reta./R + q.*cos(delta).*Reta + I2*sin(delta));
uz = -(U1/(2*pi)).*(dtilde.*q.*Reta./R + q.*sin(delta).*Reta + I4*sin(delta));
% *************************************************************************

% *** STRAINS **************************************************************
% check singularity for cos(delta)=0; % Okada (1985), equation (35) and (34)
if cos(delta) < 1E-16        % Okada (1985), equation (35)
    J1 = 0.5*alpha*(q./(R+dtilde).^2).*(2*csi.^2./(R.*(R+dtilde))-1);
    J2 = 0.5*alpha*(csi*sin(delta)./(R+dtilde).^2).*(2*q.^2./(R.*(R+dtilde))-1);
else%                       % Okada (1985), equation (34)
    J1 = alpha*(1/cos(delta))*(csi.^2./(R.*(R+dtilde).^2) - 1./(R+dtilde)) - tan(delta)*K3; 
    J2 = alpha*(1/cos(delta))*(csi.*ytilde./(R.*(R+dtilde).^2)) - tan(delta)*K1;
end;
    J3 = alpha*(-csi.*Reta./R) - J2;
    J4 = alpha*(-cos(delta)./R - q*sin(delta).*Reta./R) -J1;
    
% strike slip strains: Okada (1985), equation (31)
duxx = (U1/(2*pi)).*(csi.^2.*q.*An - J1*sin(delta));
duxy = (U1/(2*pi)).*(csi.^3.*dtilde./(R.^3.*(eta.^2+q.^2)) - (csi.^3.*An+J2)*sin(delta));
duyx = (U1/(2*pi)).*(csi.*q*cos(delta)./R.^3 + (csi.*q.^2.*An-J2)*sin(delta));
duyy = (U1/(2*pi)).*(ytilde.*q*cos(delta)./R.^3 +...
       (q.^3.*An*sin(delta) - 2*q*sin(delta).*Reta./R - (csi.^2+eta.^2)*cos(delta)./R.^3 -J4)*sin(delta));
% *************************************************************************   

% *** TILTS ***************************************************************
% strike slip tilts: Okada (1985), equation (37)
duzx = (U1/(2*pi)).*(-csi.*q.^2.*An*cos(delta) + (csi.*q./R.^3-K1)*sin(delta));
duzy = (U1/(2*pi)).*(dtilde.*q*cos(delta)./R.^3 +...
       (csi.^2.*q.*An*cos(delta) - sin(delta)./R + ytilde.*q./R.^3 -K2)*sin(delta));
% *************************************************************************   

