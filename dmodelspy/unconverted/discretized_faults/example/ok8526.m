function [ux uy uz duxx duxy duyx duyy duzx duzy] = ok8526(csi,eta,q,mu,lambda,U2,delta)
%function [ux uy uz duxx duxy duyx duyy duzx duzy ] = ok8526(csi,eta,q,mu,lambda,U2,delta)
% dip slip fault, Green's function (forward model) based on eq. (26) of
% Okada (1985). All paremeters are in SI units
% csi       dummy coordinate, defined in (22) and used in the integral (23)
% eta       dummy coordinate, defined in (22) and used in the integral (23) 
% q         projected coordinate, see (12) 
% mu        shear modulus
% lambda    Lame's first elastic parameter
% delta     dip angle from horizontal reference surface (90 degrees = vertical fault)
% U2        dip slip (U2 > 0 : the southern section of the fault goes up)
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
% Okada (1985), equation (30)
ytilde = eta*cos(delta) + q*sin(delta);
dtilde = eta*sin(delta) - q*cos(delta);
R2 = csi.^2 + ytilde.^2 + dtilde.^2; R = sqrt(R2);
X2 = csi.^2 + q.^2; X = sqrt(X2);

alpha = mu/(lambda+mu);
Rcsi = 1./(R+csi);

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

% dip slip displacement: Okada (1985), equation (26)
ux = -(U2/(2*pi)).*(q./R - I3*sin(delta)*cos(delta));
if abs(q) < 1E-16           % check singularity condition (i), pp 1148, Okada (1985)
    uy = -(U2/(2*pi)).*(ytilde.*q.*Rcsi./R + - I1*sin(delta)*cos(delta));
    uz = -(U2/(2*pi)).*(dtilde.*q.*Rcsi./R + - I5*sin(delta)*cos(delta));
else
    uy = -(U2/(2*pi)).*(ytilde.*q.*Rcsi./R + cos(delta)*atan(csi.*eta./(q.*R)) - I1*sin(delta)*cos(delta));
    uz = -(U2/(2*pi)).*(dtilde.*q.*Rcsi./R + sin(delta)*atan(csi.*eta./(q.*R)) - I5*sin(delta)*cos(delta));
end;
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
    
% dip slip strains: Okada (1985), equation (32)
duxx = (U2/(2*pi)).*(csi.*q./R.^3 + J3*sin(delta)*cos(delta));
duxy = (U2/(2*pi)).*(ytilde.*q./R.^3 - sin(delta)./R + J1*sin(delta)*cos(delta));
duyx = (U2/(2*pi)).*(ytilde.*q./R.^3 + q.*cos(delta).*Reta./R + J1*sin(delta)*cos(delta));
duyy = (U2/(2*pi)).*(ytilde.^2.*q.*Ac -...
       (2*ytilde.*Rcsi./R + csi*cos(delta).*Reta./R)*sin(delta) + J2*sin(delta)*cos(delta));
% *************************************************************************   

% *** TILTS ***************************************************************
% dip slip tilts: Okada (1985), equation (38)
duzx = (U2/(2*pi)).*(dtilde.*q./R.^3 + q.*sin(delta).*Reta./R + K3*sin(delta)*cos(delta));
duzy = (U2/(2*pi)).*(ytilde.*dtilde.*q.*Ac -...
       (2*dtilde.*Rcsi./R + csi*sin(delta).*Reta./R)*sin(delta) + K1*sin(delta)*cos(delta));
% *************************************************************************   

