% Validation of the fault model by Okada (1992)
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

clear all; close all; clc;
% define figure size
scsz = get(0,'ScreenSize');
path(pathdef);


% *************************************************************************
% *************************************************************************
% TEST 1: simple vertical fault (free surface)
faultinp = 'fault5.ss.inp';
dL = 1500; dW = 1000;
patches(faultinp,dL,dW);

% DISPLACEMENT -> Strike slip fault ***************************************
fidinp = fopen('patches.inp','r+');
M = textscan(fidinp,'%f %f %f %f %f %f %f %f %f %f %f %*[^\n]','headerlines',2);
xi = M{2}; yi = M{3}; xf = M{4}; yf = M{5}; 
Us = M{7}; Ud = M{8}; 
delta = M{9}; zt = M{10}; zb = M{11};

% Strike slip fault displacement from Coulomb 3.1
faultfid = fopen(faultinp,'r+');
M = textscan(faultfid,'%*f %f %f %f %f %*f %f %*f %f %f %f %*[^\n]','headerlines',2);
xci = M{1}; yci = M{2};
xcf = M{3}; ycf = M{4};
[xc yc zc uc vc wc] = textread('../coulomb/fault5.ss.z0.cou','%f %f %f %f %f %f','headerlines',3);
UCX = TriScatteredInterp(xc,yc,uc); UCY = TriScatteredInterp(xc,yc,vc); UCZ = TriScatteredInterp(xc,yc,wc); 

% displacement from okada85 and okada92
mu=1; nu=0.25; 
x = (0:1:25); y = x+5;   % cross-section
[u9 v9 w9] = okada92fault('strike',xi,yi,xf,yf,zt,zb,Us,delta,mu,nu,1000*xc,1000*yc,zc);
UX9 = TriScatteredInterp(1000*xc,1000*yc,u9); UY9 = TriScatteredInterp(1000*xc,1000*yc,v9); UZ9 = TriScatteredInterp(1000*xc,1000*yc,w9);

figure('Position',[100 scsz(4)/10 0.8*scsz(3) 0.5*scsz(4)],'Name','Strike slip fault (FREE SURFACE)');
subplot(1,4,[1 2]);
                plotpatches(xi,yi,xf,yf,zt,zb,delta,-7,22);
                xlabel('X (km)'); ylabel('Y (km)'); zlabel('Depth (km)');
                axis equal;
                title(sprintf('STRIKE SLIP deformation at %3.2f km',zc(1))); 
subplot(1,4,3);  
                quiver(xc,yc,u9,v9,1); hold on; 
                plot([xci xcf],[yci ycf],'r-','LineWidth',2);               % fault
                plot([0.001*xi 0.001*xf],[0.001*yi 0.001*yf],'g.');         % discretized fault
                plot(x,y,'m');                                              % cross section
                axis equal; 
                title(sprintf('Okada92: %d patches',length(xi))); 
                xlabel('X (km)'); ylabel('Y (km)');
                xlim([min(xc) max(xc)]); ylim([min(yc) max(yc)]);
subplot(1,4,4); 
                quiver(xc,yc,uc,vc,1); hold on; 
                plot([xci xcf],[yci ycf],'r-','LineWidth',2);               % fault
                plot(x,y,'m');                                              % cross section
                axis equal; title('COULOMB 3.1'); xlabel('X (km)'); ylabel('Y (km)');
                xlim([min(xc) max(xc)]); ylim([min(yc) max(yc)]);
%print -djpeg DipSlipSingleGeom
% *************************************************************************