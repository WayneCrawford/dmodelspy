% test of the 2D finite spherical source

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

% read input data for flat half-space model (McTigue, 1987)
[~, z0] = textread('inputdata.txt','%s %f',1);                              % depth
[~, P_G] = textread('inputdata.txt','%s %f',1,'headerlines',1);             % dimensionless pressure
[~, a] = textread('inputdata.txt','%s %f',1,'headerlines',2);               % radius
[~, nu] = textread('inputdata.txt','%s %f',1,'headerlines',3);              % Poisson's ratio

% load FEM model
[Xu U] = textread('FEM/udispl-0m.txt','%f %f','headerlines',1);
[Xv V] = textread('FEM/vdispl-0m.txt','%f %f','headerlines',1);
[Xw W] = textread('FEM/wdispl-0m.txt','%f %f','headerlines',1);

% DISPLACEMENT ************************************************************
% compute deformation at the surface
z = 0; x = (0.0000001:100:5000); y = x;                                % model parameters
[u3 v3 w3] = mctigue3Ddispl(0,0,z0,P_G,a,nu,x,y,z);                              % 3D numerical model
[u v w dwdx dwdy] = mctigue2D(0,0,z0,P_G,a,nu,x,y);                         % 2D analytical model

% compare 3D, 2D and FEM models
figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3]);
subplot(1,3,1); plot(x,u3,x,u,'r.',Xu,U,'g--'); xlabel('x (m)'); ylabel('x-displacement (m) at z=0');
subplot(1,3,2); plot(x,v3,x,v,'r.',Xv,V,'g--'); xlabel('x (m)'); ylabel('y-displacement (m) at z=0');
subplot(1,3,3); plot(x,w3,x,w,'r.',Xw,W,'g--'); xlabel('x (m)'); ylabel('z-displacement (m) at z=0');
legend('3D','2D','FEM'); legend('boxoff');
%filename = 'compare2D3DFEM'; print('-djpeg',filename); 


% GROUND TILT ************************************************************* 
% load FE model
[X1 dWdx] = textread('FEM/dwdx-0m.txt','%f %f','headerlines',1);
[X2 dWdy] = textread('FEM/dwdy-0m.txt','%f %f','headerlines',1);

% compare ground tilt against FE model 
figure('Position',[10 scrsz(4)/5 scrsz(3)/3 scrsz(4)/4]);
subplot(1,2,1); plot(x,dwdx,'r.',X1,dWdx,'g--'); xlabel('x (m)'); ylabel('ground tilt (E)'); legend('2D','FEM'); legend('boxoff');
subplot(1,2,2); plot(x,dwdy,'r.',X2,dWdy,'g--'); xlabel('x (m)'); ylabel('ground tilt (N)'); 
%filename = 'comparetiltFEM'; print('-djpeg',filename); 


% PRINT TABLE *************************************************************
fid = fopen('Table1.txt','w+');
M = [x'/1000 y'/1000 1000*u' 1000*v' 1000*w' 1E6*dwdx' 1E6*dwdy' 1000*u3' 1000*v3' 1000*w3'];
fprintf(fid,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',M');
