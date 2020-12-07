% Test of the 3D finite spheroidal source
% TEST3 - rotated spheroid, free surface

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


clear all; close all; clc; path(pathdef);
scrsz = get(0,'ScreenSize');

% ellipsoid parameters
x0 = 0; y0 = 0; z0 = 3000;                                                  % source location (x,y,positive depth) [m]
a = 1000;                                                                   % semi-major axis [m]  
A = 0.5;                                                                    % geometric aspect ratio
P_G = 0.1;                                                                  % excess pressure    
mu = 9.6E9;                                                                 % shear modulus
nu = 0.25;                                                                  % Poisson's ratio
theta = 0.01;                                                               % dip angle
x = (-30000:500:30000); y=x; z=zeros(size(x));                              % profile at the free surface


phi = 90;                                                                   % strike angle
[u v w] = yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z);
% plot displacement
figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3]); XLim = 10;
[X U] = textread('FEM/Rotated spheroid/u_A05_t00_p90.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,1); plot(0.001*x,u,'r.',0.001*X,U,'b'); 
title('z=0, \theta=0�, \phi=90�, A=0.5'); xlabel('X (km)'); ylabel('U (m)'); xlim([-XLim XLim]);
%legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
[X V] = textread('FEM/Rotated spheroid/v_A05_t00_p90.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,2); plot(0.001*x,v,'r.',0.001*X,V,'b'); 
title('z=0, \theta=0�, \phi=90�, A=0.5'); xlabel('X (km)'); ylabel('V (m)'); xlim([-XLim XLim]);
legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
[X W] = textread('FEM/Rotated spheroid/w_A05_t00_p90.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,3); plot(0.001*x,w,'r.',0.001*X,W,'b'); 
title('z=0, \theta=0�, \phi=90�, A=0.5'); xlabel('X (km)'); ylabel('W (m)');  xlim([-XLim XLim]);
%legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');

phi = 60;                                                                   % strike angle
[u v w] = yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z);
% plot displacement
figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3]); XLim = 10;
[X U] = textread('FEM/Rotated spheroid/u_A05_t00_p60.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,1); plot(0.001*x,u,'r.',0.001*X,U,'b'); 
title('z=0, \theta=0�, \phi=60�, A=0.5'); xlabel('X (km)'); ylabel('U (m)'); xlim([-XLim XLim]);
%legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
[X V] = textread('FEM/Rotated spheroid/v_A05_t00_p60.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,2); plot(0.001*x,v,'r.',0.001*X,V,'b'); 
title('z=0, \theta=0�, \phi=60�, A=0.5'); xlabel('X (km)'); ylabel('V (m)'); xlim([-XLim XLim]);
legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
[X W] = textread('FEM/Rotated spheroid/w_A05_t00_p60.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,3); plot(0.001*x,w,'r.',0.001*X,W,'b'); 
title('z=0, \theta=0�, \phi=60�, A=0.5'); xlabel('X (km)'); ylabel('W (m)');  xlim([-XLim XLim]);
%legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');

phi = 30;                                                                   % strike angle
[u v w] = yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z);
% plot displacement
figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3]); XLim = 10;
[X U] = textread('FEM/Rotated spheroid/u_A05_t00_p30.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,1); plot(0.001*x,u,'r.',0.001*X,U,'b'); 
title('z=0, \theta=0�, \phi=30�, A=0.50'); xlabel('X (km)'); ylabel('U (m)'); xlim([-XLim XLim]);
%legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
[X V] = textread('FEM/Rotated spheroid/v_A05_t00_p30.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,2); plot(0.001*x,v,'r.',0.001*X,V,'b'); 
title('z=0, \theta=0�, \phi=30�, A=0.50'); xlabel('X (km)'); ylabel('V (m)'); xlim([-XLim XLim]);
legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
[X W] = textread('FEM/Rotated spheroid/w_A05_t00_p30.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,3); plot(0.001*x,w,'r.',0.001*X,W,'b'); 
title('z=0, \theta=0�, \phi=30�, A=0.50'); xlabel('X (km)'); ylabel('W (m)');  xlim([-XLim XLim]);
%legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
fclose('all');
disp('That''s all folks!');

phi = 0.0;                                                                  % strike angle
[u v w] = yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z);
% plot displacement
figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3]); XLim = 10;
[X U] = textread('FEM/Rotated spheroid/u_A05_t00_p00.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,1); plot(0.001*x,u,'r.',0.001*X,U,'b'); 
title('z=0, \theta=0�, \phi=0�, A=0.5'); xlabel('X (km)'); ylabel('U (m)'); xlim([-XLim XLim]);
%legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
[X V] = textread('FEM/Rotated spheroid/v_A05_t00_p00.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,2); plot(0.001*x,v,'r.',0.001*X,V,'b'); 
title('z=0, \theta=0�, \phi=0�, A=0.5'); xlabel('X (km)'); ylabel('V (m)'); xlim([-XLim XLim]);
legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
[X W] = textread('FEM/Rotated spheroid/w_A05_t00_p00.txt','%f %f','headerlines',1);   % load FEM solution
subplot(1,3,3); plot(0.001*x,w,'r.',0.001*X,W,'b'); 
title('z=0, \theta=0�, \phi=0�, A=0.5'); xlabel('X (km)'); ylabel('W (m)');  xlim([-XLim XLim]);
%legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');


fclose('all');
disp('That''s all folks!');