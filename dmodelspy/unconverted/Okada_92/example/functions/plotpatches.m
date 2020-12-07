function plotpatches(xi,yi,xf,yf,zt,zb,delta,az,el) 
% plot a rectangle in 3D
% coordinates of corners (xi,yi,zt), (xf,yf,zt), (x0i,y0i,zb), (x0f,y0f,zb)
% xi        x start
% yi        y start
% xf        x finish
% yf        y finish
% zt        top
%           (positive downward and defined as depth below the reference surface)
% zb        bottom; zb >
%           (positive downward and defined as depth below the reference surface)
% az        view azimuth
% el        view elevation
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


    if nargin==7                                                            % set the default angles for view
        az=0; el = 0;                                                       % XZ is the default view plane 
    end;

    xi = 0.001*xi;                                                          % scale all the positions to meters
    xf = 0.001*xf;                                                          % convenient for plot
    yi = 0.001*yi;
    yf = 0.001*yf;
    zt = 0.001*zt;
    zb = 0.001*zb;
    
    delta = pi*delta/180;                                                   % angles are scaled to radians
    phi = atan((xf-xi)./(yf-yi));

    x0i = xi + (zb-zt).*cot(delta).*cos(phi);                               % lower left corner
    y0i = yi - (zb-zt).*cot(delta).*sin(phi);

    x0f = xf + (zb-zt).*cot(delta).*cos(phi);                               % lower right corner
    y0f = yf - (zb-zt).*cot(delta).*sin(phi);
    
    zm = zt+0.5*(zb-zt);                                                    % mid point
    xm = xi+0.5*(xf-xi) + zm.*cot(delta).*cos(phi);
    ym = yi - zm.*cot(delta).*sin(phi);
    
    for i=1:length(xi)                                                      % plot rectangle
        X = [xi(i) xf(i) x0f(i) x0i(i) xi(i)];                              % corners
        Y = [yi(i) yf(i) y0f(i) y0i(i) yi(i)];
        Z = [zt(i) zt(i) zb(i) zb(i) zt(i)];
        color = unifrnd(0,1,[1 3]);                                         % assign a random color
        plot3(X,Y,Z,'Color',color,'LineWidth',0.25); hold on;                            
        set(gca,'ZDir','reverse','ZLim',[0 max(zb)]);                       % reverse the Z axis
        %text(xm(i),ym(i),zm(i),num2str(i),'Color',color,'FontSize',8);      % print the rectangle number
    end;
    view(az,el);                                                            % define the 3D view