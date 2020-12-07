function [L MX MY MZ MU] = laplacian(xi,yi,xf,yf,zt,zb,U,delta,P,toggle)
% This script computes the Laplacian L of a discretized rectangular fault. 
% All variable/parameters are in SI units, angles are in degrees
% 
% OUTPUT
% L     Laplacian (P x N matrix)
%
%
% SUB-DISLOCATION PARAMETERS (all column vectors with the same length) 
% xi        x start                                     [m]
% yi        y start                                     [m]
% xf        x finish                                    [m]
% yf        y finish                                    [m]
% zt        top                                         [m]
%           (positive downward and defined as depth below the reference surface)
% zb        bottom; zb > zt                             [m]
%           (positive downward and defined as depth below the reference surface)
% delta     dip angle from horizontal reference surface (90° = vertical fault)
%           delta can be between 0° and 90° but must be different from zero!
% P         number of rows 
% % toggle  toggle [top bot] to define if the fault is breaking the surface [top = 1] or
%           a dike is opening [bot = 1]; default is [0 0];
%
% *************************************************************************
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


% ASSIGN PARAMETERS *******************************************************
% values of the top and bottom toggle
if nargin == 9
    toggle = [0 0];
end;
    top = toggle(1);
    bot = toggle(2);

NP = length(xi);                                                            % number of patches
N = round(NP/P);                                                            % number of columns
if NP == 1
    error('single rectangular fault: cannot compute smoothing operator');
end;

delta = pi*delta/180;                                                       % dip angle in radians
phi = atan2((xf-xi),(yf-yi));                                               % strike angle in radians

% compute the location of the center point of each sub-dislocation (aka sub-fault, tile or patch)
xxm = xi + 0.5*(xf-xi);                                                     % needed to compute center of subdislocation
yym = yf + (yf-yi).*(xxm-xf)./(xf-xi);                                      % needed to compute center of subdislocation
    zbm = zt + 0.5*(zb-zt);                                                 % center of subdislocation (Up)
    xfm = xxm + (zbm-zt).*cot(delta).*cos(phi);                             % center of subdislocation (East)
    yfm = yym - (zbm-zt).*cot(delta).*sin(phi);                             % center of subdislocation (North)
xfm = xfm'; yfm = yfm'; zbm = zbm';                                         % transpose vectors
% *************************************************************************


% build matrices with East, North, Up components of the center and slip of the sub-dislocation
I2 = 0;                                                                     % initialize index

for n=1:N
        I2 = I2 + P;                                                        % assign value to index                                                 
        I1 = I2 - P;                                                        % assign value to index
    MX(:,n) = xfm(1+I1:I2);                                                 % East location of center
    MY(:,n) = yfm(1+I1:I2);                                                 % North location of center
    MZ(:,n) = zbm(1+I1:I2);                                                 % Up location of center
    MU(:,n) = U(1+I1:I2);                                                   % slip
end;

% compute Laplacian L (internal subdislocation)
L = zeros(P,N);
for i=2:P-1                                                                 % loop over rows                                                             
    for j=2:N-1                                                             % loop over columns
       t4h = sqrt((MX(i,j)-MX(i-1,j))^2+(MY(i,j)-MY(i-1,j))^2+(MZ(i,j)-MZ(i-1,j))^2);
       t3h = sqrt((MX(i,j+1)-MX(i,j))^2+(MY(i,j+1)-MY(i,j))^2+(MZ(i,j+1)-MZ(i,j))^2);
       t2h = sqrt((MX(i+1,j)-MX(i,j))^2+(MY(i+1,j)-MY(i,j))^2+(MZ(i+1,j)-MZ(i,j))^2);
       t1h = sqrt((MX(i,j)-MX(i,j-1))^2+(MY(i,j)-MY(i,j-1))^2+(MZ(i,j)-MZ(i,j-1))^2);
       L(i,j) = 2*(MU(i,j-1)/(t1h*(t1h+t3h)) + MU(i+1,j)/(t2h*(t2h+t4h)) +...
                   MU(i,j+1)/(t3h*(t1h+t3h)) + MU(i-1,j)/(t4h*(t2h+t4h)) -...
                    (1/(t1h*t3h)+1/(t2h*t4h))*MU(i,j));
    end; 
end;

% top edge
i=1;                                                                        % loop over rows 
  if top == 1
    for j=1:N                                                               % loop over columns
        t2h = sqrt((MX(i+1,j)-MX(i,j))^2+(MY(i+1,j)-MY(i,j))^2+(MZ(i+1,j)-MZ(i,j))^2);
        L(i,j) = (MU(i+1,j)- MU(i,j))/t2h;
    end; 
  else
      for j=2:N-1                                                           % loop over columns
       t2h = sqrt((MX(i+1,j)-MX(i,j))^2+(MY(i+1,j)-MY(i,j))^2+(MZ(i+1,j)-MZ(i,j))^2);
       t4h = t2h;
       t3h = sqrt((MX(i,j+1)-MX(i,j))^2+(MY(i,j+1)-MY(i,j))^2+(MZ(i,j+1)-MZ(i,j))^2);
       t1h = sqrt((MX(i,j)-MX(i,j-1))^2+(MY(i,j)-MY(i,j-1))^2+(MZ(i,j)-MZ(i,j-1))^2);
       L(i,j) = 2*(MU(i,j-1)/(t1h*(t1h+t3h)) + MU(i+1,j)/(t2h*(t2h+t4h)) +...
                   MU(i,j+1)/(t3h*(t1h+t3h))                              -...
                    (1/(t1h*t3h)+1/(t2h*t4h))*MU(i,j));
      end; 
      for j=1                                                               % loop over columns
       t2h = sqrt((MX(i+1,j)-MX(i,j))^2+(MY(i+1,j)-MY(i,j))^2+(MZ(i+1,j)-MZ(i,j))^2);
       t4h = t2h;
       t3h = sqrt((MX(i,j+1)-MX(i,j))^2+(MY(i,j+1)-MY(i,j))^2+(MZ(i,j+1)-MZ(i,j))^2);
       t1h = t3h;
       L(i,j) = 2*(MU(i+1,j)/(t2h*(t2h+t4h)) +...
                   MU(i,j+1)/(t3h*(t1h+t3h))                              -...
                    (1/(t1h*t3h)+1/(t2h*t4h))*MU(i,j));
      end; 
      for j=N                                                               % loop over columns
       t2h = sqrt((MX(i+1,j)-MX(i,j))^2+(MY(i+1,j)-MY(i,j))^2+(MZ(i+1,j)-MZ(i,j))^2);
       t4h = t2h;
       t1h = sqrt((MX(i,j)-MX(i,j-1))^2+(MY(i,j)-MY(i,j-1))^2+(MZ(i,j)-MZ(i,j-1))^2);
       t3h = t1h;
       L(i,j) = 2*(MU(i,j-1)/(t1h*(t1h+t3h)) + MU(i+1,j)/(t2h*(t2h+t4h)) +...
                                                                         -...
                    (1/(t1h*t3h)+1/(t2h*t4h))*MU(i,j));
      end; 

  end;

% bottom edge
i=P;                                                                        % loop over rows  
  if bot == 1
    for j=1:N                                                               % loop over columns
       t4h = sqrt((MX(i,j)-MX(i-1,j))^2+(MY(i,j)-MY(i-1,j))^2+(MZ(i,j)-MZ(i-1,j))^2);
       L(i,j) = (MU(i,j)- MU(i-1,j))/t4h;
    end; 
  else
    for j=2:N-1                                                             % loop over columns
       t4h = sqrt((MX(i,j)-MX(i-1,j))^2+(MY(i,j)-MY(i-1,j))^2+(MZ(i,j)-MZ(i-1,j))^2);
       t3h = sqrt((MX(i,j+1)-MX(i,j))^2+(MY(i,j+1)-MY(i,j))^2+(MZ(i,j+1)-MZ(i,j))^2);
       t2h = t4h;
       t1h = sqrt((MX(i,j)-MX(i,j-1))^2+(MY(i,j)-MY(i,j-1))^2+(MZ(i,j)-MZ(i,j-1))^2);
       L(i,j) = 2*(MU(i,j-1)/(t1h*(t1h+t3h))                             +...
                   MU(i,j+1)/(t3h*(t1h+t3h)) + MU(i-1,j)/(t4h*(t2h+t4h)) -...
                    (1/(t1h*t3h)+1/(t2h*t4h))*MU(i,j));
    end;
    for j=1                                                                 % loop over columns
       t4h = sqrt((MX(i,j)-MX(i-1,j))^2+(MY(i,j)-MY(i-1,j))^2+(MZ(i,j)-MZ(i-1,j))^2);
       t3h = sqrt((MX(i,j+1)-MX(i,j))^2+(MY(i,j+1)-MY(i,j))^2+(MZ(i,j+1)-MZ(i,j))^2);
       t2h = t4h;
       t1h = t3h;
       L(i,j) = 2*(                                                      +...
                   MU(i,j+1)/(t3h*(t1h+t3h)) + MU(i-1,j)/(t4h*(t2h+t4h)) -...
                    (1/(t1h*t3h)+1/(t2h*t4h))*MU(i,j));
    end;
    for j=N                                                                 % loop over columns
       t4h = sqrt((MX(i,j)-MX(i-1,j))^2+(MY(i,j)-MY(i-1,j))^2+(MZ(i,j)-MZ(i-1,j))^2);
       t1h = sqrt((MX(i,j)-MX(i,j-1))^2+(MY(i,j)-MY(i,j-1))^2+(MZ(i,j)-MZ(i,j-1))^2);
       t3h = t1h;
       t2h = t4h;
       t1h = sqrt((MX(i,j)-MX(i,j-1))^2+(MY(i,j)-MY(i,j-1))^2+(MZ(i,j)-MZ(i,j-1))^2);
       L(i,j) = 2*(MU(i,j-1)/(t1h*(t1h+t3h))                             +...
                   MU(i-1,j)/(t4h*(t2h+t4h)) -...
                    (1/(t1h*t3h)+1/(t2h*t4h))*MU(i,j));
    end;
  end;

% left edge
for i=2:P-1                                                                 % loop over rows                                                             
    for j=1                                                                 % loop over columns
       t4h = sqrt((MX(i,j)-MX(i-1,j))^2+(MY(i,j)-MY(i-1,j))^2+(MZ(i,j)-MZ(i-1,j))^2);
       t3h = sqrt((MX(i,j+1)-MX(i,j))^2+(MY(i,j+1)-MY(i,j))^2+(MZ(i,j+1)-MZ(i,j))^2);
       t2h = sqrt((MX(i+1,j)-MX(i,j))^2+(MY(i+1,j)-MY(i,j))^2+(MZ(i+1,j)-MZ(i,j))^2);
       t1h = t3h; 
       L(i,j) = 2*(                            MU(i+1,j)/(t2h*(t2h+t4h)) +...
                   MU(i,j+1)/(t3h*(t1h+t3h)) + MU(i-1,j)/(t4h*(t2h+t4h)) -...
                    (1/(t1h*t3h)+1/(t2h*t4h))*MU(i,j));
    end; 
end;

% right edge
for i=2:P-1                                                                 % loop over rows                                                             
    for j=N                                                             % loop over columns
       t1h = sqrt((MX(i,j)-MX(i,j-1))^2+(MY(i,j)-MY(i,j-1))^2+(MZ(i,j)-MZ(i,j-1))^2);
       t2h = sqrt((MX(i+1,j)-MX(i,j))^2+(MY(i+1,j)-MY(i,j))^2+(MZ(i+1,j)-MZ(i,j))^2);
       t3h = t1h;
       t4h = sqrt((MX(i,j)-MX(i-1,j))^2+(MY(i,j)-MY(i-1,j))^2+(MZ(i,j)-MZ(i-1,j))^2);
       L(i,j) = 2*(MU(i,j-1)/(t1h*(t1h+t3h)) + MU(i+1,j)/(t2h*(t2h+t4h)) +...
                                             + MU(i-1,j)/(t4h*(t2h+t4h)) -...
                    (1/(t1h*t3h)+1/(t2h*t4h))*MU(i,j));
    end; 
end;