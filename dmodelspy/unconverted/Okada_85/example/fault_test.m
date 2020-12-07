% Validation of the fault model by Okada (1985)

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

% *************************************************************************
% Y. Okada (1985). Surface deformation due to shear and tensile faults in a
% half-space. Bull Seism Soc Am 75(4), 1135-1154.
% *************************************************************************
clear all; close all; clc;
% define figure size
scsz = get(0,'ScreenSize');
% load path
path(pathdef);

% *************************************************************************
% *************************************************************************
% TEST 1. Model validation against Table 2 from Okada (1985) **************
delta=70; x=2000; y=3000; 
x0=0; y0=0; z0=4000; mu=1; nu=0.25; U1=-1; U2=1; U3=1; phi=90; W=2000; L=3000;
xi = x0 - W*cos(pi*delta/180)*cos(pi*phi/180); xf = xi + L;
yi = y0 + W*cos(pi*delta/180)*sin(pi*phi/180); yf = yi;
zb = z0;    zt = zb - W*sin(pi*delta/180);

% print Table
fid1 = fopen('fault_test.rsl','w+');
fprintf(fid1,'Validation against Table 2 from Okada (1985)\n');
fprintf(fid1,'-----------------------------------------------------------------------------\n');
fprintf(fid1,'      \t ux \t uy \t uz \t dux/dx \t dux/dy \t duy/dx \t duy/dy \t duz/dx \t duz/dy\n');
fprintf(fid1,'-----------------------------------------------------------------------------\n');
fprintf(fid1,'Case 2\n'); 
[u v w dwdx dwdy , ~, ~, ~, duxx duxy duyx duyy] = okada85('Strike',xi,yi,xf,yf,zt,zb,U1,delta,mu,nu,x,y);
fprintf(fid1,'Strike  %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e\n',...
              u, v, w, duxx, duxy, duyx, duyy, dwdx, dwdy);
[u v w dwdx dwdy , ~, ~, ~, duxx duxy duyx duyy] = okada85('DIP',xi,yi,xf,yf,zt,zb,U2,delta,mu,nu,x,y);          
fprintf(fid1,'Dip     %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e\n',...
              u, v, w, duxx, duxy, duyx, duyy, dwdx, dwdy);
[u v w dwdx dwdy , ~, ~, ~, duxx duxy duyx duyy] = okada85('tensile',xi,yi,xf,yf,zt,zb,U3,delta,mu,nu,x,y);
fprintf(fid1,'Tensile %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e \t %5.3e\n',...
              u, v, w, duxx, duxy, duyx, duyy, dwdx, dwdy);
fprintf(fid1,'\n'); % *****************************************************         
fprintf(fid1,'Case 3\n'); delta=90; x=0; y=0;
xi = x0 - W*cot(pi*delta/180)*cos(pi*phi/180); xf = xi + L;
yi = y0 + W*cot(pi*delta/180)*sin(pi*phi/180); yf = yi;
zb = z0;    zt = zb - W*sin(pi*delta/180);
[u v w dwdx dwdy , ~, ~, ~, duxx duxy duyx duyy] = okada85('Strike',xi,yi,xf,yf,zt,zb,U1,delta,mu,nu,x,y);
fprintf(fid1,'Strike  %4.3f \t %5.3e \t %4.3f \t %4.3f \t %5.3e \t %5.3e \t %4.3f \t %4.3f \t %5.3e\n',...
              u, v, w, duxx, duxy, duyx, duyy, dwdx, dwdy);
[u v w dwdx dwdy , ~, ~, ~, duxx duxy duyx duyy] = okada85('dip',xi,yi,xf,yf,zt,zb,U2,delta,mu,nu,x,y);
fprintf(fid1,'Dip     %4.3f \t %4.3f \t %4.3f \t %4.3f \t %5.3e \t %4.3f \t %4.3f \t %4.3f \t %5.3e\n',...
              u, v, w, duxx, duxy, duyx, duyy, dwdx, dwdy);
[u v w dwdx dwdy , ~, ~, ~, duxx duxy duyx duyy] = okada85('tensile',xi,yi,xf,yf,zt,zb,U3,delta,mu,nu,x,y);
fprintf(fid1,'Tensile %5.3e \t %4.3f \t %5.3e \t %5.3e \t %4.3f \t %4.3f \t %5.3e \t %5.3e \t %4.3f\n',...
              u, v, w, duxx, duxy, duyx, duyy, dwdx, dwdy);
fprintf(fid1,'\n'); % *****************************************************
fprintf(fid1,'Case 4\n'); delta=-89.9; x=0; y=0; U=-1;
xi = x0 - W*cot(pi*delta/180)*cos(pi*phi/180); xf = xi + L;
yi = y0 + W*cot(pi*delta/180)*sin(pi*phi/180); yf = yi;
zb = z0;    zt = zb - abs(W*sin(pi*delta/180));
[u v w dwdx dwdy , ~, ~, ~, duxx duxy duyx duyy] = okada85('Strike',xi,yi,xf,yf,zt,zb,U1,delta,mu,nu,x,y);
fprintf(fid1,'Strike  %4.3f \t %5.3e \t %4.3f \t %4.3f \t %5.3e \t %5.3e \t %4.3f \t %4.3f \t %5.3e\n',...
              u, v, w, duxx, duxy, duyx, duyy, dwdx, dwdy);
[u v w dwdx dwdy , ~, ~, ~, duxx duxy duyx duyy] = okada85('dip',xi,yi,xf,yf,zt,zb,U2,delta,mu,nu,x,y);
fprintf(fid1,'Dip     %4.3f \t %4.3f \t %4.3f \t %4.3f \t %5.3e \t %4.3f \t %4.3f \t %4.3f \t %5.3e\n',...
              u, v, w, duxx, duxy, duyx, duyy, dwdx, dwdy);
[u v w dwdx dwdy , ~, ~, ~, duxx duxy duyx duyy] = okada85('tensile',xi,yi,xf,yf,zt,zb,U3,delta,mu,nu,x,y);
fprintf(fid1,'Tensile %5.3e \t %4.3f \t %5.3e \t %5.3e \t %4.3f \t %4.3f \t %5.3e \t %5.3e \t %4.3f\n',...
              u, v, w, duxx, duxy, duyx, duyy, dwdx, dwdy);
fprintf(fid1,'-----------------------------------------------------------------------------\n');
fprintf(fid1,'Case 2: Finite Source: x=2; y=3; z0=4; delta=70; L=3; W=2; mu=1; nu=0.25; U1=1; phi=90;\n'); 
fprintf(fid1,'Case 3: Finite Source: x=0; y=0; z0=4; delta=90; L=3; W=2; mu=1; nu=0.25; U1=1; phi=90;\n'); 
fprintf(fid1,'Case 4: Finite Source: x=0; y=0; z0=4; delta=-89.9; L=3; W=2; mu=1; nu=0.25; U1=1; phi=90;\n'); 

fclose(fid1); 
% *************************************************************************
% *************************************************************************



% *************************************************************************
% *************************************************************************
% TEST 2: Validation agaisnt COULOMB 3.1 
% NOTE: in Coulomb the fault geometry is in km!!!!!!!

% DISPLACEMENT -> Strike slip fault ***************************************
% strike-slip fault displacement from Coulomb 3.1
[xc yc zc uc vc wc] = textread('coulomb/okada85ss.cou','%f %f %f %f %f %f','headerlines',3);
UCX = TriScatteredInterp(xc,yc,uc); UCY = TriScatteredInterp(xc,yc,vc); UCZ = TriScatteredInterp(xc,yc,wc); 

% displacement
xi = 4000; yi = 4000; xf = 25000; yf = 20000; delta=60; mu=1; nu=0.25; U=1; zt =1000; zb = 10000;
[u v w dwdx dwdy]  = okada85('Strike',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,1000*xc,1000*yc);
UX = TriScatteredInterp(1000*xc,1000*yc,u); UY = TriScatteredInterp(1000*xc,1000*yc,v); UZ = TriScatteredInterp(1000*xc,1000*yc,w);
DWDX = TriScatteredInterp(1000*xc,1000*yc,dwdx); DWDY = TriScatteredInterp(1000*xc,1000*yc,dwdy);

figure('Position',[100 scsz(4)/10 0.5*scsz(3) 0.5*scsz(4)],'Name','Strike slip fault');
subplot(1,2,1);  
                quiver(xc,yc,u,v,1); hold on; 
                plot([0.001*xi 0.001*xf],[0.001*yi 0.001*yf],'r-','LineWidth',2);
                plot([-50 50],[0 0],'g-','LineWidth',1);
                axis equal; title('Okada85 - strike slip'); xlabel('X, in kilometers'); ylabel('Y, in kilometers');
                xlim([min(xc) max(xc)]); ylim([min(yc) max(yc)])
subplot(1,2,2); 
                quiver(xc,yc,uc,vc,1); hold on; 
                plot([0.001*xi 0.001*xf],[0.001*yi 0.001*yf],'r-','LineWidth',2);
                plot([-50 50],[0 0],'g-','LineWidth',1);
                axis equal; title('COULOMB - strike slip'); xlabel('X, in kilometers'); ylabel('Y, in kilometers');
                xlim([min(xc) max(xc)]); ylim([min(yc) max(yc)])

figure('Position',[100 scsz(4)/10 0.66*scsz(3) 0.33*scsz(4)],'Name','Strike slip fault');              
x = (-50:1:50); y = zeros(size(x));
subplot(1,3,1);  
                plot(x,UX(1000*x,1000*y),'.',x,UCX(x,y),'r'); 
                ylabel('u, in meters - East'); xlabel('x, in kilometers'); 
subplot(1,3,2); 
                plot(x,UY(1000*x,1000*y),'.',x,UCY(x,y),'r'); legend('boxoff'); legend('Okada','Coulomb'); legend('boxoff'); legend('Location','Best');
                ylabel('v, in meters - North'); xlabel('x, in kilometers'); 
subplot(1,3,3); 
                plot(x,UZ(1000*x,1000*y),'.',x,UCZ(x,y),'r'); 
                ylabel('w, in meters - Up'); xlabel('x, in kilometers'); 
% *************************************************************************
%[x' UX(1000*x,1000*y)' UY(1000*x,1000*y)' UZ(1000*x,1000*y)' 1E5*DWDX(1000*x,1000*y)' 1E5*DWDY(1000*x,1000*y)'];
% *************************************************************************


% DISPLACEMENT -> Dip slip fault ***************************************
% dip-slip fault displacement from Coulomb 3.1
[xc yc zc uc vc wc] = textread('coulomb/okada85ds.cou','%f %f %f %f %f %f','headerlines',3);
UCX = TriScatteredInterp(xc,yc,uc); UCY = TriScatteredInterp(xc,yc,vc); UCZ = TriScatteredInterp(xc,yc,wc); 

% displacement
xi = 4000; yi = 4000; xf = 25000; yf = 20000; delta=60; mu=1; nu=0.25; U=1; zt =1000; zb = 10000;
[u v w dwdx dwdy]  = okada85('dip',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,1000*xc,1000*yc);
UX = TriScatteredInterp(1000*xc,1000*yc,u); UY = TriScatteredInterp(1000*xc,1000*yc,v); UZ = TriScatteredInterp(1000*xc,1000*yc,w);
DWDX = TriScatteredInterp(1000*xc,1000*yc,dwdx); DWDY = TriScatteredInterp(1000*xc,1000*yc,dwdy);

figure('Position',[100 scsz(4)/10 0.5*scsz(3) 0.5*scsz(4)],'Name','Dip slip fault');
subplot(1,2,1);  
                quiver(xc,yc,u,v,1); hold on; 
                plot([0.001*xi 0.001*xf],[0.001*yi 0.001*yf],'r-','LineWidth',2);
                plot([-50 50],[0 0],'g-','LineWidth',1);
                axis equal; title('Okada85 - dip slip'); xlabel('X, in kilometers'); ylabel('Y, in kilometers');
                xlim([min(xc) max(xc)]); ylim([min(yc) max(yc)])
subplot(1,2,2); 
                quiver(xc,yc,uc,vc,1); hold on; 
                plot([0.001*xi 0.001*xf],[0.001*yi 0.001*yf],'r-','LineWidth',2);
                plot([-50 50],[0 0],'g-','LineWidth',1);
                axis equal; title('COULOMB - dip slip'); xlabel('X, in kilometers'); ylabel('Y, in kilometers');
                xlim([min(xc) max(xc)]); ylim([min(yc) max(yc)])

figure('Position',[100 scsz(4)/10 0.66*scsz(3) 0.33*scsz(4)],'Name','Dip slip fault');              
x = (-50:1:50); y = zeros(size(x));
subplot(1,3,1);  
                plot(x,UX(1000*x,1000*y),'.',x,UCX(x,y),'r'); 
                ylabel('u, in meters - East'); xlabel('x, in kilometers'); 
subplot(1,3,2); 
                plot(x,UY(1000*x,1000*y),'.',x,UCY(x,y),'r'); 
                ylabel('v, in meters - North'); xlabel('x, in kilometers'); 
subplot(1,3,3); 
                plot(x,UZ(1000*x,1000*y),'.',x,UCZ(x,y),'r'); 
                ylabel('w, in meters - Up'); xlabel('x, in kilometers'); legend('boxoff'); legend('Okada','Coulomb'); legend('Location','Best');
% *************************************************************************
%[x' UX(1000*x,1000*y)' UY(1000*x,1000*y)' UZ(1000*x,1000*y)' 1E5*DWDX(1000*x,1000*y)' 1E5*DWDY(1000*x,1000*y)'];
% *************************************************************************


% DISPLACEMENT -> Tensile fault ***************************************
clear UCX UCY UCZ xc yc;

% Tensile fault displacement from FEM
[xfe yfe ufe] = textread('FEM/tensile/u_xy.txt','%f %f %f','headerlines',1); UCX = TriScatteredInterp(xfe,yfe,ufe); 
[xfe yfe vfe] = textread('FEM/tensile/v_xy.txt','%f %f %f','headerlines',1); UCY = TriScatteredInterp(xfe,yfe,vfe); 

% displacement from Okada85
xi = -5000; yi = 0; xf = 5000; yf = 0; delta=90; mu=1; nu=0.25; U=1; zt =1000; zb = 6000;
x = (2*xi:1000:2*xf); y = x; [xc,yc] = meshgrid(x,y);
[u v w]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xc,yc);

figure('Position',[100 scsz(4)/10 0.5*scsz(3) 0.5*scsz(4)],'Name','Tensile fault');
subplot(1,2,1);  
                quiver(0.001*xc,0.001*yc,u,v,1); hold on; 
                plot([-10 10],[-10 10],'g-','LineWidth',1);
                plot([0.001*xi 0.001*xf],[0.001*yi 0.001*yf],'r-','LineWidth',2);
                axis equal; title('Okada85 - tensile opening'); xlabel('X, in kilometers'); ylabel('Y, in kilometers');
subplot(1,2,2);  
                quiver(0.001*xc,0.001*yc,UCX(xc,yc),UCY(xc,yc),1); hold on; 
                plot([-10 10],[-10 10],'g-','LineWidth',1);
                plot([0.001*xi 0.001*xf],[0.001*yi 0.001*yf],'r-','LineWidth',2);
                axis equal; title('FEM - tensile opening'); xlabel('X, in kilometers'); ylabel('Y, in kilometers');


figure('Position',[100 scsz(4)/10 0.66*scsz(3) 0.33*scsz(4)],'Name','Tensile fault');              
clear xfe ufe vfe wfe y u v w

subplot(1,3,1);  
    [xfe ufe] = textread('FEM/tensile/u_y.txt','%f %f','headerlines',1); 
    [u , ~, ~]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xfe,xfe);
    plot(0.001*xfe,u,'.',0.001*xfe,ufe,'r'); 
    ylabel('u, in meters - East'); xlabel('x, in kilometers'); 
subplot(1,3,2); clear xfe;
    [xfe vfe] = textread('FEM/tensile/v_y.txt','%f %f','headerlines',1); 
    [~, v, ~]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xfe,xfe);
    plot(0.001*xfe,v,'.',0.001*xfe,vfe,'r'); 
    ylabel('v, in meters - North'); xlabel('x, in kilometers'); 
subplot(1,3,3); clear xfe;
    [xfe wfe] = textread('FEM/tensile/w_y.txt','%f %f','headerlines',1); 
    [~, ~, w]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xfe,xfe);
    plot(0.001*xfe,w,'.',0.001*xfe,wfe,'r'); 
    ylabel('w, in meters - Up'); xlabel('x, in kilometers'); 
    legend('Okada','FEM'); legend('boxoff'); legend('Location','Best');
% *************************************************************************
[~, ~, ~, dwdx dwdy]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xfe,xfe);
[1E-3*xfe 1E-3*xfe u v w 1E5*dwdx 1E5*dwdy]
% *************************************************************************


% TILT -> Tensile fault ***************************************
figure('Position',[100 scsz(4)/10 0.5*scsz(3) 0.5*scsz(4)],'Name','Tensile fault');
clear xfe 

subplot(1,2,1);  
    [xfe dwdxfe] = textread('FEM/tensile/dwdx_y.txt','%f %f','headerlines',1); 
    [~, ~, ~, dwdx dwdy]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xfe,xfe);
    plot(0.001*xfe,dwdx,'b',0.001*xfe,dwdxfe,'r.'); 
    title('TILT: \tau_E = dw/dx'); ylabel('EAST component'); xlabel('x, in kilometers'); 
subplot(1,2,2); clear xfe;
    [xfe dwdyfe] = textread('FEM/tensile/dwdy_y.txt','%f %f','headerlines',1); 
     [~, ~, ~, dwdx dwdy]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xfe,xfe);
    plot(0.001*xfe,dwdy,'b',0.001*xfe,dwdyfe,'r.'); 
    title('TILT: \tau_N = dw/dy'); ylabel('NORTH component'); xlabel('x, in kilometers'); 
    legend('Okada','FEM'); legend('boxoff'); legend('Location','Best');
% *************************************************************************
% *************************************************************************


% STRAIN -> Tensile fault ***************************************
figure('Position',[100 scsz(4)/10 0.66*scsz(3) 0.33*scsz(4)],'Name','STRAIN Tensile fault');              
clear xfe ufe vfe wfe y u v w

subplot(1,3,1);  
    [xfe eeafe] = textread('FEM/tensile/eea_y.txt','%f %f','headerlines',1); 
    [~, ~, ~, ~, ~, eea]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xfe,xfe);
    plot(0.001*xfe,eea,'b',0.001*xfe,eeafe,'r.'); 
    ylabel('e_a'); xlabel('y, in kilometers'); 
subplot(1,3,2); clear xfe;
    [xfe gamma1fe] = textread('FEM/tensile/gamma1_y.txt','%f %f','headerlines',1); 
    [~, ~, ~, ~, ~, ~, gamma1]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xfe,xfe);
    plot(0.001*xfe,gamma1,'b',0.001*xfe,gamma1fe,'r.'); 
    ylabel('\gamma_1'); xlabel('y, in kilometers'); 
subplot(1,3,3); clear xfe;
    [xfe gamma2fe] = textread('FEM/tensile/gamma2_y.txt','%f %f','headerlines',1); 
    [~, ~, ~, ~, ~, eea gamma1 gamma2]  = okada85('tensile',xi,yi,xf,yf,zt,zb,U,delta,mu,nu,xfe,xfe);
    plot(0.001*xfe,gamma2,'b',0.001*xfe,gamma2fe,'r.'); 
    ylabel('\gamma_1'); xlabel('y, in kilometers'); 
    legend('Okada','FEM'); legend('boxoff'); legend('Location','Best');
% *************************************************************************
% *************************************************************************
