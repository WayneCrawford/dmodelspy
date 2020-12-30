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

% plot displacement at z=0;
[u v w dV] = fialko01(x0,y0,z0,P_G,a,nu,x,y,0);
figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3]); XLim = 10;
[rr ur] = textread('FEM/sill_ur_r.txt','%f %f','headerlines',1);         % read FEM model for radial deformation
[rw W] = textread('FEM/sill_w_r.txt','%f %f','headerlines',1);          % read FEM model for vertical deformation
UR = spline(rr,ur,r);

subplot(1,3,1); plot(0.001*r,u,'r.',0.001*r,UR.*x./r,'b'); 
title('SILL-LIKE SOURCE'); xlabel('radial distance, in kilometers'); ylabel('U, in meters');  xlim([0 XLim]);
subplot(1,3,2); plot(0.001*r,v,'r.',0.001*r,UR.*y./r,'b'); 
title('z=0, a=1000 m, z_0=1000 m'); xlabel('radial distance, in kilometers'); ylabel('V, in meters'); xlim([0 XLim]);
subplot(1,3,3); plot(0.001*r,w,'r.',0.001*rw,W,'b'); 
xlabel('radial distance, in kilometers'); ylabel('W, in meters');  xlim([0 XLim]);
legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');


% plot displacement at z=500 m;
[u v w dV] = fialko01(x0,y0,z0,P_G,a,nu,x,y,500);
figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3]); XLim = 15;
[rr ur] = textread('FEM/sill_ur_r_500m.txt','%f %f','headerlines',1);         % read FEM model for radial deformation
[rw W] = textread('FEM/sill_w_r_500m.txt','%f %f','headerlines',1);          % read FEM model for vertical deformation
UR = spline(rr,ur,r);

subplot(1,3,1); plot(0.001*r,u,'r.',0.001*r,UR.*x./r,'b'); 
title('SILL-LIKE SOURCE'); xlabel('radial distance, in kilometers'); ylabel('U, in meters');  
xlim([0 XLim]); ylim([0 Inf]);
subplot(1,3,2); plot(0.001*r,v,'r.',0.001*r,UR.*y./r,'b'); 
title('z=500 m, a=1000 m, z_0=1000 m'); xlabel('radial distance, in kilometers'); ylabel('V, in meters'); 
xlim([0 XLim]); ylim([0 Inf]);
subplot(1,3,3); plot(0.001*r,w,'r.',0.001*rw,W,'b'); 
xlabel('radial distance, in kilometers'); ylabel('W, in meters');  xlim([0 XLim]); ylim([0 Inf]);
legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');


fclose('all');
disp('That''s all folks!');