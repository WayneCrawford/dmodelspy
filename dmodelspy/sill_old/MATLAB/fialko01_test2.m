% Test of the 3D sill-like source

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
scrsz = get(0,'ScreenSize');

% penny-crack parameters
x0 = 0; y0 = 0; z0 = 1000; 
P_G = 0.01; a = 1000; nu=0.25; 

% profile
x = 1E3*(eps:.05:5); y = 2*x; r = sqrt(x.^2+y.^2);

% TILT ********************************************************************
[~, ~, ~, ~, dwdx dwdy ] = fialko01(x0,y0,z0,P_G,a,nu,x,y,0);

% plot ground tilt
figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3]); XLim = 10;
[X DWDX] = textread('FEM/sill_dwdx_x.txt','%f %f','headerlines',1);         % read FEM model for radial deformation
[Y DWDY] = textread('FEM/sill_dwdy_y.txt','%f %f','headerlines',1);          % read FEM model for vertical deformation

subplot(1,3,1); plot(0.001*x,dwdx,'r.',0.001*X,DWDX,'b'); 
title('SILL-LIKE SOURCE'); xlabel('X, in kilometers'); ylabel('Tilt (North)'); xlim([0 XLim]);
subplot(1,3,2); plot(0.001*y,dwdy,'r.',0.001*Y,DWDY,'b'); 
title('z=0, a=1000 m, z_0=1000 m'); xlabel('Y, in kilometers'); ylabel('Tilt (East)'); xlim([0 XLim]);
legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');


% STRAIN ******************************************************************
[~, ~, ~, ~, ~, ~, eea gamma1 gamma2] = fialko01(x0,y0,z0,P_G,a,nu,x,y,500);

% plot strain
figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3]); XLim = 10;
[R EEA] = textread('FEM/sill_eea_r_500m.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,1); plot(0.001*r,eea,'r.',0.001*R,EEA,'b'); 
title('SILL-LIKE SOURCE'); xlabel('Radial distance, in kilometers'); ylabel('e_a'); xlim([0 XLim]);
[R GAMMA1] = textread('FEM/sill_gamma1_r_500m.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,2); plot(0.001*r,gamma1,'r.',0.001*R,GAMMA1,'b'); 
title('z=500 m, a=1000 m, z_0=1000 m'); xlabel('Radial distance, in kilometers'); ylabel('\gamma_1'); xlim([0 XLim]);
[R GAMMA2] = textread('FEM/sill_gamma2_r_500m.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,3); plot(0.001*r,gamma2,'r.',0.001*R,GAMMA2,'b'); 
xlabel('Radial distance, in kilometers'); ylabel('\gamma_2');  xlim([0 XLim]);
legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');


fclose('all');
disp('That''s all folks!');