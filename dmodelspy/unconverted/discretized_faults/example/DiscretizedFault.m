% This script shows how to build a forward model of a discretized fault
% using the "okada85fault" function
% The script makes use of seven functions 
%   (a) okada85fault        compute the ground displacement, tilt and strain
%       (i) okada85         implement the Okada (1985) rectangular dislocation model
%           (1) ok8525      implement the strike slip fault model (Okada,1985)
%           (2) ok8526      implement the dip slip fault model (Okada,1985)    
%           (3) ok8527      implement the tensile dislocation model (Okada,1985)
%   (b) patches             discretized a segment in rectangular tiles (patches)
%       (i) segments        divide a fault in segments

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
faultinp = 'fault5.ss.inp';                                                 % read the input file with the fault geometry
dL = 1000; dW = 1500;                                                       % set the grid size for length (dL) and width (dW) 
patches(faultinp,dL,dW);                                                    % call function

% STEP 2: input the geometry of the discretized fault 
fidinp = fopen('patches.inp','r+');                                         % read the input file with the discretized fault geometry
M = textscan(fidinp,'%f %f %f %f %f %f %f %f %f %f %f %*[^\n]','headerlines',2);
xi = M{2}; yi = M{3}; xf = M{4}; yf = M{5}; Us = M{7}; Ud = M{8}; delta = M{9}; zt = M{10}; zb = M{11};

% STEP 3: compute displacement
[xc yc] = meshgrid((-25E3:5E3:50E3),(-25E3:5E3:50E3));                      % define benchmarks location
mu=1; nu=0.25;                                                              % define elastic parameters    
[u v w] = okada85fault('strike',xi,yi,xf,yf,zt,zb,Us,delta,mu,nu,xc,yc);    % call function


% *************************************************************************
% PLOT A FIGURE - only for display purposes *******************************

faultfid = fopen(faultinp,'r+');                                            % read the original fault geometry
M = textscan(faultfid,'%*f %f %f %f %f %*f %f %*f %f %f %f %*[^\n]','headerlines',2);
xci = 1E3*M{1}; yci = 1E3*M{2}; xcf = 1E3*M{3}; ycf = 1E3*M{4};

figure('Position',[100 scsz(4)/10 0.8*scsz(3) 0.5*scsz(4)],'Name','Strike slip fault');
subplot(1,2,1);                                                             % plot the discretized fault
                plotpatches(xi,yi,xf,yf,zt,zb,delta,-7,22);
                xlabel('X (km)'); ylabel('Y (km)'); zlabel('Depth (km)');
                title(sprintf('STRIKE SLIP: discretized fault (%d patches)',length(xi))); 
subplot(1,2,2);  
                quiver(xc,yc,u,v,0.5); hold on;                               % plot the orizontal displacement field
                plot([xci xcf],[yci ycf],'r-','LineWidth',2);
                plot([0.001*xi 0.001*xf],[0.001*yi 0.001*yf],'g.');
                axis equal; title('Okada85 (discretized fault)'); xlabel('X (km)'); ylabel('Y (km)');
                xlim([min(min(xc)) max(max(xc))]); ylim([min(min(yc)) max(max(yc))]);
% *************************************************************************