% This script validates the function okada85D.m and okadagf.m
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

% STEP 1: build a discretized fault using the function "patches"
faultinp = 'fault5.ss.inp';                                                 % read the input file with the fault geometry
dL = [2000 500 3000 1000 1500]; dW = 500;                                                      % set the grid size for length (dL) and width (dW) 
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

% plot D, Laplacian and their difference
figure('Position',[100 scsz(4)/10 0.98*scsz(3) 0.25*scsz(4)],'Name','D and L');
subplot(1,4,1); 
spy(D); axis equal; xlim([1 length(D)]); ylim([1 length(D)]);
title('Smoothing operator, D');

subplot(1,4,2); 
[X,Y] = meshgrid(0.001*xi(1:P:end),0.001*yi(1:N:end));
xlim([min(min(X)) max(max(X))]); ylim([min(min(Y)) max(max(Y))]);
pcolor(X,Y,U); colorbar; axis equal; set(gca,'YDir','reverse')
title('Slip, U'); xlabel('X (km)'); ylabel('Depth (km)');
subplot(1,4,3); 
[X,Y] = meshgrid(0.001*xi(1:P:end),0.001*yi(1:N:end));
xlim([min(min(X)) max(max(X))]); ylim([min(min(Y)) max(max(Y))]);
pcolor(X,Y,ML); colorbar; axis equal; set(gca,'YDir','reverse')
title('Numerical Laplacian, L'); xlabel('X (km)'); ylabel('Depth (km)');
subplot(1,4,4); 
[X,Y] = meshgrid(0.001*xi(1:P:end),0.001*yi(1:N:end));
xlim([min(min(X)) max(max(X))]); ylim([min(min(Y)) max(max(Y))]);
pcolor(X,Y,ML-L); colorbar; axis equal; set(gca,'ZDir','reverse')
title('L-D*U'); xlabel('X (km)'); ylabel('Depth (km)');



% compute a numerical example
clear all;
faultinp = 'fault3.ss.inp';                                                 % read the input file with the fault geometry
dL = [1000 500 3000]; dW = 5000;                                                      % set the grid size for length (dL) and width (dW) 
[xi yi xf yf Us Ud delta zt zb P] = patches(faultinp,dL,dW);                % call function

toggle = [0 0];
D  = okada85D(xi,yi,xf,yf,zt,zb,delta,P,toggle);                            % call function

Ltmp = D*Ud'; 
Utmp = Ud';
N=length(xi)/P;
for n=1:N
    L(:,n) = Ltmp(1+(n-1)*P:n*P);
    U(:,n) = Utmp(1+(n-1)*P:n*P);
end;

[ML MX MY MZ MU] = laplacian(xi,yi,xf,yf,zt,zb,Ud,delta,P,toggle);

figure
spy(D,'sr',8); 
%pcolor(D); 
%image(D);
%colorbar; 
axis equal; set(gca,'YDir','reverse')
axis equal; xlim([1 length(D)]); ylim([1 length(D)]);
title('Smoothing operator, D');


% compute weight with Jessica's matrix (good only for a homogenous grid)
%[Lap, Lap_inv] = Modelwt(P, N, dL, dW,0);
%LJtmp = Lap*Us';
%for n=1:N
%    LJ(:,n) = LJtmp(1+(n-1)*P:n*P);
%end;

% *************************************************************************
% *************************************************************************

% STEP 3: compute Green's function
xc = [500 3000 5000]; yc = [1500 10000 8000];                               % define benchmarks location
mu=1; nu=0.25;                                                              % define elastic parameters    
[ugf vgf wgf] = okada85gf('strike',xi,yi,xf,yf,zt,zb,delta,mu,nu,xc,yc);    % call function

us = ugf*Ud';                                                               % compute displacement
vs = vgf*Ud';
ws = wgf*Ud';

% STEP 4: compute displacement (forward model)
mu=1; nu=0.25;                                                              % define elastic parameters    
[u v w] = okada85fault('strike',xi,yi,xf,yf,zt,zb,Ud,delta,mu,nu,xc,yc);    % call function

% compare displacements
[u' us v' vs w' ws]

% *************************************************************************
% *************************************************************************


% *************************************************************************
% Plot a figure of the grid - only for display purposes *******************
figure;
plotpatches(xi,yi,xf,yf,zt,zb,delta,-7,22); hold on;
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Depth (km)');
title(sprintf('Discretized fault (%d patches)',length(xi))); 


