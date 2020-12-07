% COMPUTE A NUMERICAL EXAMPLE
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


clear all; close all; clc;
% define figure size
scsz = get(0,'ScreenSize');
% load path
path(pathdef);

% STEP 1: build a discretized fault using the function "patches"
faultinp = 'fault3.ss.inp';                                                 % read the input file with the fault geometry
dL = [3000 1000 3000]; dW = 5000;                                                      % set the grid size for length (dL) and width (dW) 
[xi yi xf yf Us Ud delta zt zb P] = patches(faultinp,dL,dW);                % call function

% STEP 2: compute smoothing operator D ************************************
toggle = [0 0];
D  = okada85D(xi,yi,xf,yf,zt,zb,delta,P,toggle);                            % call function

% compute Laplacian D*U
Ltmp = D*Us'; 
Utmp = Us';
N=length(xi)/P;
for n=1:N
    L(:,n) = Ltmp(1+(n-1)*P:n*P);
    U(:,n) = Utmp(1+(n-1)*P:n*P);
end;

% compute the laplacian in a completely different way for testing
[ML MX MY MZ MU] = laplacian(xi,yi,xf,yf,zt,zb,Us,delta,P,toggle);

% Plot a figure of the grid - only for display purposes *******************
plotpatches(xi,yi,xf,yf,zt,zb,delta,-7,22); hold on;
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Depth (km)');
title(sprintf('Discretized fault (%d tiles)',length(xi))); 

